"""Summarize repeated seeds of each ordered query/receiver and preserve all three paired comparisons."""

from collections import defaultdict
import csv
from pathlib import Path
from statistics import mean
import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent


def main():
    groups = defaultdict(list)
    for panel, path in [
        ("fixed", HERE / "results/cells.csv"),
        ("BRD4", BENCH / "search_trial_nov28/results/cells.csv"),
        ("4OCZ", BENCH / "search_trial_affinity/controls/four_versions/cells.csv"),
    ]:
        with path.open() as f:
            for r in csv.DictReader(f):
                groups[
                    panel, r["cohort"], r["target"], r["query"], r["receiver"], r["arm"]
                ].append(r)
    output = []
    for key, rows in sorted(groups.items()):
        valid = [r for r in rows if r["status"] == "ok"]
        row = dict(
            zip(("panel", "cohort", "target", "query", "receiver", "arm"), key),
            planned=len(rows),
            valid=len(valid),
            empty=sum(r["status"] == "empty" for r in rows),
            failed=sum(r["status"] not in ("ok", "empty") for r in rows),
        )
        first = [float(r["first"]) for r in valid]
        if first:
            quant = np.quantile(
                first, [0.25, 0.5, 0.75, 0.90, 0.95], method="averaged_inverted_cdf"
            )
            row.update(zip(("q1", "median", "q3", "p90", "p95"), quant))
            for metric in ("first", "best5", "best", "regret"):
                row["mean_" + metric] = mean(float(r[metric]) for r in valid)
        for arm in ("official", "parent", "stable"):
            ref = {
                r["seed"]: float(r["first"])
                for r in groups[(*key[:-1], arm)]
                if r["status"] == "ok"
            }
            diff = [
                float(r["first"]) - ref[r["seed"]] for r in valid if r["seed"] in ref
            ]
            row["paired_" + arm] = len(diff)
            row["gain_" + arm] = sum(v < -1e-8 for v in diff)
            row["tie_" + arm] = sum(abs(v) <= 1e-8 for v in diff)
            row["loss_" + arm] = sum(v > 1e-8 for v in diff)
        output.append(row)
    with (HERE / "results/by_direction.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            list(dict.fromkeys(k for row in output for k in row)),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()
