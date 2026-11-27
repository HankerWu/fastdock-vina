"""Pair the new regression with the frozen results per seed and crystal copy."""

import argparse
import csv
import re
from pathlib import Path

import numpy as np


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def compare_runs(old_run, new_run, out):
    before = {
        (r["method"], r["receptor"], r["compound"], r["seed"]): r
        for r in read(old_run / "runs.csv")
    }
    after = {
        (r["method"], r["receptor"], r["compound"], r["seed"]): r
        for r in read(new_run / "runs.csv")
    }
    assert before.keys() == after.keys() and len(before) == 72
    rows = []
    for key, prev in before.items():
        next_row = after[key]
        log = "_".join(key) + ".log"
        old_log = (old_run / "pools" / log).read_text()
        new_log = (new_run / "pools" / log).read_text()
        old_search = [
            l for l in old_log.splitlines() if l.startswith("Search trajectory:")
        ]
        new_search = [
            l for l in new_log.splitlines() if l.startswith("Search trajectory:")
        ]
        row = dict(zip(("method", "receptor", "compound", "seed"), key))
        row.update(
            old_score=prev["score"],
            new_score=next_row["score"],
            old_modes=prev["modes"],
            new_modes=next_row["modes"],
            same_pool_bytes=prev["pool_sha256"] == next_row["pool_sha256"],
            timing_source="frozen_control" if key[0] == "official" else "new_run",
            old_wall=prev["wall_s"],
            new_wall=next_row["wall_s"],
            old_cpu=prev["cpu_s"],
            new_cpu=next_row["cpu_s"],
            mc_counters_equal=old_search == new_search if old_search else "",
        )
        for stage in ("hunt_evals", "second_evals", "explicit_evals"):
            for prefix, text in (("old", old_log), ("new", new_log)):
                counts = re.findall(stage + r"=(\d+)", text)
                row[prefix + "_" + stage] = sum(map(int, counts)) if counts else ""
        rows.append(row)
    write(out / "solver_runs.csv", rows)


def compare(old, new, summary):
    keys = ("compound", "reference", "chain", "primary", "receptor", "seed", "method")
    get_key = lambda r: tuple(r[k] for k in keys)
    before = {get_key(r): r for r in read(old / "geometry_cells.csv")}
    after = {get_key(r): r for r in read(new / "geometry_cells.csv")}
    assert before.keys() == after.keys()
    pairs = []
    for k, base in before.items():
        out = dict(zip(keys, k))
        for metric in ("top1", "best5", "oracle", "regret"):
            lo, hi = (float(row[metric]) for row in (base, after[k]))
            out.update(
                {f"old_{metric}": lo, f"new_{metric}": hi, f"delta_{metric}": hi - lo}
            )
        pairs.append(out)
    write(summary / "solver_pairs.csv", pairs)
    stats = []
    for method in ("official", "pro", "double"):
        rows = [r for r in pairs if r["primary"] == "True" and r["method"] == method]
        assert len(rows) == 16
        for metric in ("top1", "best5", "oracle", "regret"):
            deltas = [r["delta_" + metric] for r in rows]
            q25, med, q75, p90, p95 = np.percentile(deltas, [25, 50, 75, 90, 95])
            stats.append(
                {
                    "method": method,
                    "metric": metric,
                    "count": len(rows),
                    "old_mean": np.mean([r["old_" + metric] for r in rows]),
                    "new_mean": np.mean([r["new_" + metric] for r in rows]),
                    "mean_delta": np.mean(deltas),
                    "q25": q25,
                    "median_delta": med,
                    "q75": q75,
                    "p90": p90,
                    "p95": p95,
                    "improved": sum(x < -1e-12 for x in deltas),
                    "worse": sum(x > 1e-12 for x in deltas),
                    "equal": sum(abs(x) <= 1e-12 for x in deltas),
                }
            )
    write(summary / "solver_contrasts.csv", stats)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("old", "new", "old_run", "new_run", "out"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    compare(a.old, a.new, a.out)
    compare_runs(a.old_run, a.new_run, a.out)
