"""Aggregate each ordered receptor pair across its repeated seeds, keeping all four raw endpoints."""

import collections
import csv
import sys
from pathlib import Path
from statistics import mean

import numpy as np

path = Path(sys.argv[1])
with (path / "cells.csv").open() as stream:
    rows = list(csv.DictReader(stream))
groups = collections.defaultdict(list)
for r in rows:
    key = tuple(r[k] for k in ("cohort", "target", "query", "receiver", "arm"))
    groups[key].append(r)
results = []
for key, group in sorted(groups.items()):
    valid = [r for r in group if r["status"] == "ok"]
    first = [float(r["first"]) for r in valid]
    row = dict(zip(("cohort", "target", "query", "receiver", "arm"), key))
    row.update(
        planned=len(group),
        valid=len(valid),
        empty=sum(r["status"] == "empty" for r in group),
        failed=sum(r["status"] not in ("ok", "empty") for r in group),
    )
    q = (
        np.quantile(first, [0.25, 0.5, 0.75, 0.9, 0.95], method="averaged_inverted_cdf")
        if first
        else [None] * 5
    )
    row.update(zip(("q1", "median", "q3", "p90", "p95"), q))
    for metric in ("first", "best5", "best", "regret"):
        row["mean_" + metric] = mean(float(r[metric]) for r in valid) if valid else None
    for arm in ("official", "parent", "stable"):
        paired = {
            r["seed"]: float(r["first"])
            for r in groups[(*key[:-1], arm)]
            if r["status"] == "ok"
        }
        diffs = [
            float(r["first"]) - paired[r["seed"]] for r in valid if r["seed"] in paired
        ]
        row["pairs_" + arm] = len(diffs)
        row["gain_" + arm] = sum(d < -1e-8 for d in diffs)
        row["tie_" + arm] = sum(abs(d) <= 1e-8 for d in diffs)
        row["loss_" + arm] = sum(d > 1e-8 for d in diffs)
    results.append(row)
with (path / "by_direction.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=results[0].keys(), lineterminator="\n")
    writer.writeheader()
    writer.writerows(results)
