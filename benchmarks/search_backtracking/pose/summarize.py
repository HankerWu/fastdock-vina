"""Write per-target, cap, ligand and receiver pose outcomes and paired changes."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

from backtracking_panel import ARMS, paired


def summarize(rows):
    def hits(col):
        return sum(float(r[col]) <= 2 for r in rows)

    def total(col):
        return sum(int(r[col]) for r in rows)

    out = {
        "cells": len(rows),
        "top1": hits("top1_A"),
        "top5": hits("best5_A"),
        "returned_oracle": hits("oracle_A"),
        "pool_misses": len(rows) - hits("oracle_A"),
        "ranking_misses": hits("oracle_A") - hits("top1_A"),
        "median_top1_A": median(float(r["top1_A"]) for r in rows),
        "top1_with_clashes": sum(int(r["top1_clash_pairs"]) != 0 for r in rows),
        "seconds": sum(float(r["seconds"]) for r in rows),
        "joint_evals": sum(
            total(k) for k in ("hunt_evals", "second_evals", "explicit_evals")
        ),
        "local_calls": sum(
            total(k) for k in ("hunt_calls", "second_calls", "explicit_calls")
        ),
    }
    for k in (
        "hunt_evals",
        "second_evals",
        "explicit_evals",
        "hunt_rejected",
        "hunt_exhausted",
        "save_rejected",
        "save_exhausted",
        "hunt_searches",
        "save_searches",
    ):
        out[k] = total(k) if rows[0][k] != "" else ""
    return out


def run(work):
    with (work / "cells.csv").open() as f:
        rows = list(csv.DictReader(f))
    group = defaultdict(list)
    for r in rows:
        for dim, label in (
            ("all", "all"),
            ("target", r["target"]),
            ("target_cap", r["target"] + "/" + r["cap"]),
            ("cap", r["cap"]),
            ("kind", r["kind"]),
            ("case", r["query"] + "/" + r["receiver"]),
        ):
            group[dim, label, r["arm"]].append(r)
    paired.write(
        work / "summary.csv",
        [
            dict(group=dim, label=label, arm=arm, **summarize(v))
            for (dim, label, arm), v in sorted(group.items())
        ],
    )
    cells = defaultdict(dict)
    for r in rows:
        cells[r["key"].rsplit("_", 1)[0]][r["arm"]] = r
    changes = []
    for key, arms in sorted(cells.items()):
        assert set(arms) == set(ARMS)
        for ref, test in (
            ("parent", "halving"),
            ("parent", "quadratic"),
            ("halving", "quadratic"),
        ):
            a, b = arms[ref], arms[test]

            def delta(col, a=a, b=b):
                return float(b[col]) - float(a[col])

            def hitdelta(col, a=a, b=b):
                return int(float(b[col]) <= 2) - int(float(a[col]) <= 2)

            changes.append(
                {
                    "cell": key,
                    "target": a["target"],
                    "ref": ref,
                    "candidate": test,
                    "top1_delta_A": delta("top1_A"),
                    "top1_delta_hit": hitdelta("top1_A"),
                    "top5_delta_hit": hitdelta("best5_A"),
                    "oracle_delta_hit": hitdelta("oracle_A"),
                    "delta_joint_evals": int(
                        sum(
                            delta(k)
                            for k in ("hunt_evals", "second_evals", "explicit_evals")
                        )
                    ),
                    "delta_seconds": delta("seconds"),
                }
            )
    paired.write(work / "contrasts.csv", changes)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    run(p.parse_args().work)
