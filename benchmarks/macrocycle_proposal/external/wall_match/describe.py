"""Summarize the preselected caps on all validation seeds, with actual process costs."""

import argparse
import csv
from pathlib import Path
from statistics import mean

from matched_time import SEEDS, TARGETS, write


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main(path):
    runs = rows(path / "all_runs.csv")
    geometry = {
        (r["pdb"], int(r["seed"]), r["arm"]): r for r in rows(path / "structures.csv")
    }
    table = []
    for target in TARGETS:
        for seed in SEEDS:
            rr = {
                (int(r["rep"]), r["arm"]): r
                for r in runs
                if r["pdb"] == target and int(r["seed"]) == seed
            }
            parent = [int(rr[i, "parent"]["elapsed_ns"]) / 1e6 for i in range(4)]
            lock = [int(rr[i, "lock"]["elapsed_ns"]) / 1e6 for i in range(4)]
            gp, gl = [geometry[target, seed, a] for a in ("parent", "lock")]
            er = int(gl["joint"]) / int(gp["joint"])
            wr = mean(lock) / mean(parent)
            stats = {
                "pdb": target,
                "seed": seed,
                "parent_cap": gp["cap"],
                "lock_cap": gl["cap"],
                "parent_wall_ms": round(mean(parent), 4),
                "lock_wall_ms": round(mean(lock), 4),
                "wall_ratio": round(wr, 6),
                "parent_joint": gp["joint"],
                "lock_joint": gl["joint"],
                "eval_ratio": round(er, 6),
                "both_costs_within_parent": er <= 1 and wr <= 1,
                "within_95_percent_wall": wr <= 0.95,
                "parent_top1_A": gp["top1_A"],
                "lock_top1_A": gl["top1_A"],
                "parent_oracle_A": gp["oracle_A"],
                "lock_oracle_A": gl["oracle_A"],
                "parent_top1_cg_A": gp["top1_cg_distance_A"],
                "lock_top1_cg_A": gl["top1_cg_distance_A"],
                "parent_top1_clashes": gp["top1_clashes"],
                "lock_top1_clashes": gl["top1_clashes"],
            }
            stats.update(
                {f"wall_ratio_r{i}": round(lock[i] / parent[i], 6) for i in range(4)}
            )
            table.append(stats)
    write(path / "paired.csv", table)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("validation", type=Path)
    a = ap.parse_args()
    main(a.validation)
