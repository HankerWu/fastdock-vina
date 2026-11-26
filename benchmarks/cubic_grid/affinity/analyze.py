"""Read every mode pool, rescore its top using the original topology, and compare exact Kd pairs."""

import argparse
import csv
import itertools
import re
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from run import SEEDS, STATES, box_args


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def number(text, prefix):
    line = next(l for l in text.splitlines() if l.startswith(prefix))
    return float(line.split(":", 1)[1].strip().split()[0])


def analyze(work, binary, jobs):
    prep = work / "prepared"
    panel = list(csv.DictReader((work / "panel.csv").open()))

    def read(case):
        r, seed, arm = case
        state = STATES[r["target"]]
        stem = f"{state}_{r['cid']}_{seed}_{arm}"
        text = (work / "pools" / (stem + ".pdbqt")).read_text()
        log = (work / "pools" / (stem + ".log")).read_text()
        top = text.split("ENDMDL", 1)[0]
        score = number(top, "REMARK VINA RESULT:")
        with tempfile.TemporaryDirectory() as tmp:
            ligand = Path(tmp) / "top.pdbqt"
            ligand.write_text(top.split("\n", 1)[1])
            cmd = [
                binary,
                "--receptor",
                prep / (state + ".pdbqt"),
                "--ligand",
                prep / (r["cid"] + ".pdbqt"),
                "--pose_coordinates",
                ligand,
                "--score_only",
                "--cpu",
                "1",
                "--unbound_energy",
                str(number(top, "REMARK UNBOUND:")),
                *box_args(prep, state),
            ]
            replay = subprocess.check_output(cmd, text=True)
            readback = number(replay, "Estimated Free Energy of Binding")
        counters = [
            dict(re.findall(r"(\w+)=(\d+)", l))
            for l in log.splitlines()
            if l.startswith("Search trajectory:")
        ]
        assert len(counters) == 2
        assert text.endswith("ENDMDL\n")
        return dict(
            target=r["target"],
            cid=r["cid"],
            state=state,
            seed=seed,
            arm=arm,
            score=score,
            readback=readback,
            discrepancy=round(abs(readback - score), 5),
            modes=text.count("ENDMDL"),
            wall_s=float(log.rsplit("wall_s ", 1)[1]),
            evaluations=sum(
                int(t[k]) for t in counters for k in ("hunt_evals", "second_evals")
            ),
        )

    with ThreadPoolExecutor(jobs) as ex:
        cells = list(
            ex.map(
                read,
                [(r, k, a) for r in panel for k in SEEDS for a in ("native", "cubic")],
            )
        )
    write(work / "cells.csv", cells)
    scores = defaultdict(dict)
    for c in cells:
        scores[c["target"], c["cid"], c["arm"]][c["seed"]] = c["score"]
    pairs = []
    summary = []
    for target in STATES:
        lig = [r for r in panel if r["target"] == target]
        n = len(lig)
        for k in (*SEEDS, "mean"):
            selected = SEEDS if k == "mean" else (k,)
            get = lambda r, a: np.mean(
                [scores[target, r["cid"], a][s] for s in selected]
            )
            rows = []
            for i, j in itertools.combinations(range(n), 2):
                strong, weak = sorted((i, j), key=lambda x: float(lig[x]["kd_nm"]))
                s, w = lig[strong], lig[weak]
                if float(w["kd_nm"]) < 3 * float(s["kd_nm"]):
                    continue
                row = dict(target=target, seed=k, strong=s["cid"], weak=w["cid"])
                for a in ("native", "cubic", "heavy"):
                    delta = (
                        float(s["heavy"]) - float(w["heavy"])
                        if a == "heavy"
                        else get(w, a) - get(s, a)
                    )
                    row[a] = 1 if delta > 1e-8 else 0 if delta < -1e-8 else 0.5
                rows.append(row)
            pairs.extend(rows)
            result = dict(
                target=target,
                seed=k,
                ligands=n,
                scaffolds=len({r["scaffold"] for r in lig}),
                pairs=len(rows),
            )
            for a in ("native", "cubic", "heavy"):
                result[a] = sum(row[a] for row in rows)
            result.update(
                gained=sum(r["cubic"] > r["native"] for r in rows),
                lost=sum(r["cubic"] < r["native"] for r in rows),
            )
            for a in ("native", "cubic"):
                result[a + "_rho"] = round(
                    float(
                        spearmanr(
                            [get(r, a) for r in lig], [float(r["kd_nm"]) for r in lig]
                        ).statistic
                    ),
                    4,
                )
            if k == "mean":
                rng = np.random.default_rng(650053)
                idx = {r["cid"]: i for i, r in enumerate(lig)}
                counts = rng.multinomial(n, [1 / n] * n, size=10000)
                weights = np.array(
                    [
                        counts[:, idx[r["strong"]]] * counts[:, idx[r["weak"]]]
                        for r in rows
                    ]
                )
                denominator = weights.sum(axis=0)
                num = np.array([r["cubic"] - r["native"] for r in rows]) @ weights
                diffs = num[denominator > 0] / denominator[denominator > 0]
                ci = np.quantile(diffs, [0.025, 0.975])
                result.update(delta_lo=round(ci[0], 4), delta_hi=round(ci[1], 4))
            else:
                result.update(delta_lo="", delta_hi="")
            summary.append(result)
    write(work / "pairs.csv", pairs)
    write(work / "summary.csv", summary)
    print(
        (work / "summary.csv").read_text(),
        f"Max absolute readback discrepancy {max(c['discrepancy'] for c in cells)}",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", type=int, default=2)
    a = p.parse_args()
    analyze(a.work.resolve(), a.binary.resolve(), a.jobs)
