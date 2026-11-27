"""Re-dock all eight Holt ligands with original seeds, budget and assay aggregation."""

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import itertools
from pathlib import Path
import random
import re
import statistics
import subprocess
import time
from scipy.stats import kendalltau, spearmanr
from run_panel import write


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "out", "parent", "release"):
        p.add_argument(name, type=Path)
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--confirm", action="store_true")
    a = p.parse_args()
    pool = a.out / "pools"
    pool.mkdir(parents=True, exist_ok=True)
    labels = list(
        csv.DictReader(
            (
                Path(__file__).resolve().parents[1]
                / "cubic_grid/human_panel/affinity/ligands.csv"
            ).open()
        )
    )
    recs = ("1FKG", "1FKH")
    seeds = (524309, 654323, 786433) if a.confirm else (131071, 196613, 262147)
    arms = ("parent", "release")
    jobs = list(itertools.product([l["compound"] for l in labels], recs, seeds, arms))
    random.Random(831337).shuffle(jobs)

    def run(job):
        ligand, rec, seed, arm = job
        dst = pool / f"{ligand}-{rec}-{seed}-{arm}"
        args = [
            getattr(a, arm),
            "--ligand",
            a.work / "affinity" / (ligand + ".pdbqt"),
            "--maps",
            a.work / "maps" / rec / "receptor",
            "--scoring",
            "ad4",
            "--cpu",
            "1",
            "--exhaustiveness",
            "3",
            "--max_evals",
            "24000",
            "--seed",
            str(seed),
            "--num_modes",
            "9",
            "--energy_range",
            "100",
            "--search_counters",
            "--out",
            dst.with_suffix(".pdbqt"),
        ]
        t = time.perf_counter()
        with dst.with_suffix(".log").open("w") as f:
            subprocess.run(args, stdout=f, stderr=subprocess.STDOUT, check=True)
        wall = time.perf_counter() - t
        scores = [
            float(l.split()[3])
            for l in dst.with_suffix(".pdbqt").read_text().splitlines()
            if l.startswith("REMARK VINA RESULT:")
        ]
        counts = re.findall(
            r"hunt_evals=(\d+).*?second_evals=(\d+)",
            dst.with_suffix(".log").read_text(),
        )
        assert len(counts) == 3 and len(scores) > 0
        return dict(
            ligand=ligand,
            receiver=rec,
            seed=seed,
            arm=arm,
            score=scores[0],
            modes=len(scores),
            wall_s=wall,
            evals=sum(int(t) for pair in counts for t in pair),
        )

    with ThreadPoolExecutor(a.jobs) as ex:
        rows = list(ex.map(run, jobs))
    write(a.out / "affinity_search.csv", rows)
    pred = {
        (l["compound"], r, arm): statistics.median(
            x["score"]
            for x in rows
            if (x["ligand"], x["receiver"], x["arm"]) == (l["compound"], r, arm)
        )
        for l, r, arm in itertools.product(labels, recs, arms)
    }
    for l, arm in itertools.product(labels, arms):
        pred[l["compound"], "mean", arm] = statistics.mean(
            pred[l["compound"], r, arm] for r in recs
        )
    s = []
    pairs = []
    for r, arm in itertools.product((*recs, "mean"), arms):
        points = [l for l in labels if l["ki_low_nm"] == l["ki_high_nm"]]
        measured = [float(l["ki_low_nm"]) for l in points]
        guessed = [pred[l["compound"], r, arm] for l in points]
        correct = 0
        count = 0
        for first, second in itertools.permutations(labels, 2):
            if float(second["ki_low_nm"]) < 3 * float(first["ki_high_nm"]):
                continue
            d = pred[second["compound"], r, arm] - pred[first["compound"], r, arm]
            credit = int(d > 1e-10) + 0.5 * int(abs(d) <= 1e-10)
            pairs.append(
                dict(
                    receiver=r,
                    arm=arm,
                    strong=first["compound"],
                    weak=second["compound"],
                    delta=d,
                    credit=credit,
                )
            )
            correct += credit
            count += 1
        s.append(
            dict(
                receiver=r,
                arm=arm,
                correct=correct,
                count=count,
                rho=spearmanr(measured, guessed).statistic,
                tau=kendalltau(measured, guessed).statistic,
            )
        )
    write(a.out / "affinity_search_stats.csv", s)
    write(a.out / "affinity_pairs.csv", pairs)


if __name__ == "__main__":
    main()
