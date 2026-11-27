"""Redock both defaults on the unchanged Holt FKBP12 AD4 assay panel."""

import argparse
import csv
import itertools
import json
import resource
import statistics
import sys
import tarfile
from pathlib import Path

from scipy.stats import kendalltau, spearmanr


def read(path):
    with path.open() as inp:
        return list(csv.DictReader(inp))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("repo", "work", "parent", "fixed"):
        p.add_argument(key, type=lambda s: Path(s).resolve())
    p.add_argument("--search-counters", action="store_true")
    a = p.parse_args()
    benchmark = a.repo / "benchmarks/cubic_grid"
    source = benchmark / "fkbp_ad4"
    sys.path.insert(0, str(source))
    from run_official import COMPOUNDS, RECEPTORS, SEEDS, run_one, sha, table

    work = a.work
    work.mkdir(parents=True, exist_ok=True)
    with tarfile.open(benchmark / "human_panel/affinity-pools.tar.xz") as tar:
        for i in COMPOUNDS:
            tar.extract(f"affinity/{i}.pdbqt", work, filter="data")
    meta = json.loads((source / "results/metadata.json").read_text())
    for n in COMPOUNDS:
        assert sha(work / f"affinity/{n}.pdbqt") == meta["ligands"][str(n)]
    for r in RECEPTORS:
        with tarfile.open(source / f"maps-{r}.tar.xz") as tar:
            tar.extractall(work, filter="data")
        for name, digest in meta["maps"][r].items():
            assert sha(work / "maps" / r / name) == digest
    for arm in ("parent", "fixed"):
        root = work / arm
        (root / "pools").mkdir(exist_ok=True, parents=True)
        if not (root / "maps").exists():
            (root / "maps").symlink_to(work / "maps", target_is_directory=True)

    def child_cpu():
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        return usage.ru_utime + usage.ru_stime

    rows = []
    for idx, (rec, compound, seed) in enumerate(
        itertools.product(RECEPTORS, COMPOUNDS, SEEDS)
    ):
        arms = ("parent", "fixed") if idx % 2 == 0 else ("fixed", "parent")
        for arm in arms:
            start = child_cpu()
            row = run_one(
                work,
                work / "affinity",
                work / arm,
                getattr(a, arm),
                (rec, compound, seed, "ad4"),
                search_counters=a.search_counters,
            )
            rows.append({**row, "arm": arm, "cpu_s": child_cpu() - start})
    table(work / "runs.csv", rows)
    for row in rows:
        name = f"{row['compound']}_{row['receiver']}_{row['seed']}_ad4.pdbqt"
        path = work / row["arm"] / "pools" / name
        result_lines = [
            line
            for line in path.read_text().splitlines()
            if line.startswith("REMARK VINA RESULT:")
        ]
        assert sha(path) == row["sha256"] and len(result_lines) == row["modes"]
        assert float(result_lines[0].split()[3]) == row["score"]
    labels = {
        int(r["compound"]): (float(r["ki_low_nm"]), float(r["ki_high_nm"]))
        for r in read(benchmark / "human_panel/affinity/ligands.csv")
    }
    official = [
        {**r, "arm": "official"}
        for r in read(source / "results/cells.csv")
        if r["arm"] == "ad4"
    ]
    points = [n for n, (lo, hi) in labels.items() if lo == hi]
    assert len(official) == 48 and len(rows) == 96 and len(points) == 7
    scores, pairs, stats = [], [], []
    for arm in ("official", "parent", "fixed"):
        medians = {}
        for rec, n in itertools.product(RECEPTORS, COMPOUNDS):
            subset = [
                r
                for r in [*official, *rows]
                if r["arm"] == arm and r["receiver"] == rec and int(r["compound"]) == n
            ]
            assert sorted(int(r["seed"]) for r in subset) == sorted(SEEDS)
            medians[rec, n] = statistics.median(float(r["score"]) for r in subset)
        for rec in (*RECEPTORS, "receptor_mean"):
            sc = {
                n: medians[rec, n]
                if rec in RECEPTORS
                else statistics.mean(medians[r, n] for r in RECEPTORS)
                for n in COMPOUNDS
            }
            for n, value in sc.items():
                scores.append(
                    {
                        "arm": arm,
                        "receiver": rec,
                        "compound": n,
                        "predicted": value,
                        "low_nM": labels[n][0],
                        "high_nM": labels[n][1],
                    }
                )
            order = []
            for i, j in itertools.permutations(labels, 2):
                if labels[j][0] >= 3 * labels[i][1]:
                    margin = sc[j] - sc[i]
                    credit = (margin > 0) + 0.5 * (margin == 0)
                    order.append(credit)
                    pairs.append(
                        {
                            "arm": arm,
                            "receiver": rec,
                            "stronger": i,
                            "weaker": j,
                            "margin": margin,
                            "correct": credit,
                        }
                    )
            assert len(order) == 15
            obs, pred = [labels[n][0] for n in points], [sc[n] for n in points]
            stats.append(
                {
                    "arm": arm,
                    "receiver": rec,
                    "correct": sum(order),
                    "total": len(order),
                    "points": len(points),
                    "spearman": spearmanr(obs, pred).statistic,
                    "kendall_b": kendalltau(obs, pred).statistic,
                    "epimer_8_10": sc[10] - sc[8],
                    "epimer_9_11": sc[11] - sc[9],
                }
            )
    table(work / "scores.csv", scores)
    table(work / "pairs.csv", pairs)
    table(work / "assay.csv", stats)
    costs = [
        {
            "arm": arm,
            "wall_s": sum(r["seconds"] for r in rows if r["arm"] == arm),
            "cpu_s": sum(r["cpu_s"] for r in rows if r["arm"] == arm),
            "modes": sum(r["modes"] for r in rows if r["arm"] == arm),
            "binary_sha256": sha(getattr(a, arm)),
        }
        for arm in ("parent", "fixed")
    ]
    table(work / "costs.csv", costs)


if __name__ == "__main__":
    main()
