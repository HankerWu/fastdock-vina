"""Aggregate archived runs, treating seeds and budgets as repeated conditions."""

import argparse
import csv
from pathlib import Path
import statistics


def read(p):
    with p.open() as fp:
        return list(csv.DictReader(fp))


def summarize(results):
    rows = []
    for name in ("panel", "transfer", "catd", "lower", "fresh"):
        runs = read(results / name / "structures.csv")
        timings = read(results / name / "timings.csv")
        for pdb, arm in sorted({(r["pdb"], r["arm"]) for r in runs}):
            subset = [r for r in runs if (r["pdb"], r["arm"]) == (pdb, arm)]
            keys = {r["key"] for r in subset}
            times = [t for t in timings if t["key"] in keys]
            top = sum(float(r["top1_A"]) <= 2 for r in subset)
            oracle = sum(float(r["oracle_A"]) <= 2 for r in subset)
            rows.append(
                dict(
                    experiment=name,
                    pdb=pdb,
                    arm=arm,
                    conditions=len(subset),
                    top1=top,
                    top5=sum(float(r["best5_A"]) <= 2 for r in subset),
                    pool=oracle,
                    rank1_given_hit=(top / oracle if oracle else ""),
                    mean_top1_A=statistics.mean(float(r["top1_A"]) for r in subset),
                    mean_best_A=statistics.mean(float(r["oracle_A"]) for r in subset),
                    aggregate_joint_evals=sum(int(r["joint"]) for r in subset),
                    mean_wall_ms=statistics.mean(
                        int(t["wall_ns"]) * 1e-6 for t in times
                    ),
                    mean_cpu_ms=statistics.mean(
                        float(t["cpu_s"]) * 1000 for t in times
                    ),
                    peak_KiB=max(int(t["max_rss_KiB"]) for t in times),
                )
            )
    dest = results / "summary.csv"
    with dest.open("w") as f:
        w = csv.DictWriter(f, rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results", type=Path)
    summarize(p.parse_args().results)
