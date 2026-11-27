"""Measure fixed-pool feasibility and aggregate deterministic paired runs."""

import argparse
import csv
from collections import defaultdict
from itertools import groupby
from pathlib import Path
from statistics import mean

from external_run import write


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def audit(directory, repeat):
    modes = defaultdict(list)
    for r in rows(directory / "models.csv"):
        modes[r["key"]].append(r)
    controls = []
    runs = rows(directory / "structures.csv")
    repeats = {r["key"]: r for r in rows(repeat / "runs.csv")}
    for r in runs:
        t = repeats[r["key"]]
        metrics = (
            "sha256",
            "modes",
            "joint",
            "hunt_calls",
            "hunt_evals",
            "second_calls",
            "second_evals",
            "explicit_calls",
            "explicit_evals",
        )
        if any(t[k] != r[k] for k in metrics):
            raise ValueError("Deterministic repeat failed for " + r["key"])
        r["mean_wall_ms"] = (int(r["elapsed_ns"]) + int(t["elapsed_ns"])) / 2e6
        pool = modes[r["key"]]
        valid = [m for m in pool if float(m["cg_distance_A"]) <= 3]
        top = pool[0]
        chosen = valid[0] if valid else top
        controls.append(
            {
                "key": r["key"],
                "pdb": r["pdb"],
                "arm": r["arm"],
                "seed": r["seed"],
                "cap": r["cap"],
                "size": len(pool),
                "count_under_3A": len(valid),
                "all_over_3A": int(not valid),
                "changed": int(chosen != top),
                "selected_rank": chosen["rank"],
                "before_rmsd_A": top["rmsd_A"],
                "after_rmsd_A": chosen["rmsd_A"],
                "before_cg_A": top["cg_distance_A"],
                "after_cg_A": chosen["cg_distance_A"],
            }
        )
    write(directory / "fixed_pools.csv", controls)

    def key(r):
        return r["pdb"], r["cap"], r["arm"]

    summary = []
    for (pdb, cap, arm), group in groupby(sorted(runs, key=key), key=key):
        group = list(group)
        summary.append(
            {
                "pdb": pdb,
                "cap": cap,
                "arm": arm,
                "n": len(group),
                "top1": sum(float(r["top1_A"]) <= 2 for r in group),
                "top5": sum(float(r["best5_A"]) <= 2 for r in group),
                "oracle": sum(float(r["oracle_A"]) <= 2 for r in group),
                "mean_top1_A": round(mean(float(r["top1_A"]) for r in group), 5),
                "mean_oracle_A": round(mean(float(r["oracle_A"]) for r in group), 5),
                "top1_bond_over_3A": sum(
                    float(r["top1_cg_distance_A"]) > 3 for r in group
                ),
                "mean_joint_evals": round(mean(int(r["joint"]) for r in group), 2),
                "mean_wall_ms_two_repeats": round(
                    mean(r["mean_wall_ms"] for r in group), 2
                ),
            }
        )
    write(directory / "summary.csv", summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("repeat", type=Path)
    a = p.parse_args()
    audit(a.directory, a.repeat)
