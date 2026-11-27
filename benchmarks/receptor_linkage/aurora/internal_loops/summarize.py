"""Read all arm-level outputs, retaining the eight paired conditions and full mode distributions."""

import argparse
import csv
from pathlib import Path

import numpy as np


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def save(path, rows):
    with path.open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def collect(folders, filename):
    return [
        row for p in folders if (p / filename).exists() for row in read(p / filename)
    ]


def analyze(folders, dest):
    dest.mkdir(exist_ok=True, parents=True)
    cells, modes, fixed = [
        collect(folders, name + ".csv") for name in ("cells", "modes", "fixed")
    ]
    key = lambda r: (r["arm"], r["drug"], r["seed"])
    if len({key(r) for r in cells}) != len(cells):
        raise ValueError("Overlapping docking runs")
    dry = {(r["drug"], r["seed"]): r for r in cells if r["arm"] == "dry"}
    for r in cells:
        other = dry[r["drug"], r["seed"]]
        for name in ("top1", "best5", "oracle", "regret"):
            r[name + "_minus_dry"] = float(r[name]) - float(other[name])
    for r in fixed:
        base = next(
            q
            for q in fixed
            if q["arm"] == "dry" and q["drug"] == r["drug"] and q["pose"] == r["pose"]
        )
        r["delta_score"] = round(float(r["score"]) - float(base["score"]), 3)
        r["delta_decoy_minus_witness"] = round(
            float(r["decoy_minus_witness"]) - float(base["decoy_minus_witness"]), 3
        )
    records = []
    for arm in sorted({r["arm"] for r in cells}):
        for drug in ("both", "at9283", "danusertib"):
            rows = [
                r
                for r in cells
                if r["arm"] == arm and (drug == "both" or drug == r["drug"])
            ]
            for metric in ("top1", "best5", "oracle", "regret", "top1_minus_dry"):
                vals = np.array([float(r[metric]) for r in rows])
                q = np.percentile(vals, [25, 50, 75, 90, 95])
                records.append(
                    {
                        "arm": arm,
                        "drug": drug,
                        "metric": metric,
                        "n": len(rows),
                        "mean": vals.mean(),
                        "q25": q[0],
                        "median": q[1],
                        "q75": q[2],
                        "p90": q[3],
                        "p95": q[4],
                        "maximum": vals.max(),
                        "above_10": int((vals > 10).sum()) if metric == "top1" else "",
                        "modes": sum(int(r["n"]) for r in rows),
                        "wall_s": sum(float(r["wall_s"]) for r in rows),
                        "cpu_s": sum(float(r["cpu_s"]) for r in rows),
                        "evals": sum(int(r["evals"]) for r in rows)
                        if all(r["evals"] for r in rows)
                        else "",
                    }
                )
    for name, data in (
        ("cells", cells),
        ("modes", modes),
        ("fixed", fixed),
        ("distributions", records),
    ):
        save(dest / (name + ".csv"), data)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("inputs", nargs="+", type=Path)
    args = p.parse_args()
    analyze(args.inputs, args.out)
