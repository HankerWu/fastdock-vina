"""Paired per-cell and per-target summaries from raw per-mode measurements."""

import argparse
import csv
import itertools
from collections import defaultdict
from operator import itemgetter
from pathlib import Path
from statistics import median


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def summarize(work, cap):
    lookup = defaultdict(dict)
    for arm in ("native", "cubic"):
        with (work / f"modes_{arm}_{cap}.csv").open() as f:
            modes = list(csv.DictReader(f))
        for key, g in itertools.groupby(
            modes,
            lambda r: (r["target"], r["query"], r["receiver"], r["kind"], r["seed"]),
        ):
            rows = list(g)
            lookup[key][arm] = {
                "top1": float(rows[0]["rmsd"]),
                "top5": min(float(x["rmsd"]) for x in rows[:5]),
                "oracle": min(float(x["rmsd"]) for x in rows),
                "n_modes": len(rows),
            }
    cells = []
    for key, pair in sorted(lookup.items()):
        d = dict(zip(("target", "query", "receiver", "kind", "seed"), key))
        for arm, metrics in pair.items():
            d.update({arm + "_" + k: v for k, v in metrics.items()})
        cells.append(d)
    write(work / "cells.csv", cells)
    groups = {"all": cells}
    for key in ("kind", "target"):
        for name, g in itertools.groupby(
            sorted(cells, key=itemgetter(key)), itemgetter(key)
        ):
            groups[name] = list(g)
    result = []
    for name, group in groups.items():
        for metric in ("top1", "top5", "oracle"):
            arr = [(r["native_" + metric], r["cubic_" + metric]) for r in group]
            native = sum(a <= 2 for a, b in arr)
            cubic = sum(b <= 2 for a, b in arr)
            gains = sum(a > 2 and b <= 2 for a, b in arr)
            losses = sum(a <= 2 and b > 2 for a, b in arr)
            result.append(
                {
                    "subset": name,
                    "metric": metric,
                    "n_cells": len(group),
                    "native_hits": native,
                    "cubic_hits": cubic,
                    "gains": gains,
                    "losses": losses,
                    "native_median_A": median(a for a, b in arr),
                    "cubic_median_A": median(b for a, b in arr),
                    "median_paired_change_A": median(b - a for a, b in arr),
                }
            )
    write(work / "summary.csv", result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--cap", type=int, default=24000)
    a = p.parse_args()
    summarize(a.work, a.cap)
