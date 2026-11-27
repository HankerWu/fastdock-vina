"""Tabulate dependent panel cells, grouped by family, cap, label geometry and direction."""

import argparse
import csv
import statistics
from pathlib import Path
from inputs import reference


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def report(work):
    rows = read(work / "cells.csv")
    timing = read(work / "timing_replays.csv")
    costs = {}
    for row in timing:
        costs.setdefault(row["key"], []).append(row)
    groups = {}
    clear = {
        (r["query"], r["receiver"]): int(r["reference_clashes"]) == 0
        for r in read(work / "reference_geometry.csv")
    }
    for r in rows:
        group = [
            "all",
            r["target"],
            r["kind"],
            "cap_" + r["cap"],
            "label_clear" if clear[r["query"], r["receiver"]] else "label_clash",
            r["query"] + "_to_" + r["receiver"],
        ]
        for label in group:
            groups.setdefault((label, r["arm"]), []).append(r)
    totals = []
    for (label, arm), cells in sorted(groups.items()):
        row = dict(group=label, arm=arm, cells=len(cells))
        for col, k in [
            ("top1_A", "top1_hits"),
            ("best5_A", "top5_hits"),
            ("oracle_A", "returned_hits"),
        ]:
            values = [float(c[col]) for c in cells]
            row[k] = sum(v <= 2 for v in values)
            row["median_" + col] = statistics.median(values)
        row["conditional_rank1_success"] = (
            row["top1_hits"] / row["returned_hits"] if row["returned_hits"] else ""
        )
        for k in (
            "hunt_calls",
            "hunt_evals",
            "second_calls",
            "second_evals",
            "explicit_inputs",
            "explicit_calls",
            "explicit_evals",
            "modes",
            "top1_clash_pairs",
            "skipped",
            "checks",
            "records",
            "eligible",
            "replaced",
        ):
            row[k] = sum(int(c[k]) for c in cells) if cells[0][k] else ""
        row["joint_evals"] = sum(
            row[k] for k in ("hunt_evals", "second_evals", "explicit_evals")
        )
        row["seconds"] = sum(float(c["seconds"]) for c in cells)
        for k in ("seconds", "cpu_seconds"):
            row["sum_per_cell_median_" + k] = (
                sum(
                    statistics.median(float(t[k]) for t in costs[c["key"]])
                    for c in cells
                )
                if arm != "upstream"
                else ""
            )
        totals.append(row)
    reference.write(work / "summary.csv", totals)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    report(p.parse_args().work)
