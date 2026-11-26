"""Allocate FastDock's measured CPU saving to search, independently of RMSD."""
import argparse
import csv
import json
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import replay

CAL_SEED = 74317
PER_COMPLEX_SEEDS = (93529, 93563)
LOW, HIGH = 12000, 96000


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("protocol", choices=("cognate", "reciprocal"))
    for arg in ("fixtures", "output", "official", "branch", "analytic"):
        p.add_argument(arg, type=Path)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--per-complex-calibration", type=Path,
                   help="Use previously recorded calibration.csv CPU times for individual complex budgets and fresh seeds")
    args = p.parse_args()
    fixtures, output = args.fixtures.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    binary = {name: exe.resolve() for name, exe in [
        ("official", args.official), ("native", args.branch), ("geometry_analytic", args.analytic)]}
    sources = list(replay.experiments(args.protocol, fixtures))
    calibration = []
    with ProcessPoolExecutor(args.jobs) as pool:
        if args.per_complex_calibration:
            with args.per_complex_calibration.open() as stream:
                calibration = list(csv.DictReader(stream))
            seeds = PER_COMPLEX_SEEDS
        else:
            for arm, cap in (("official", LOW), ("native", LOW), ("native", HIGH)):
                raw = output / "calibration" / (arm + str(cap))
                raw.mkdir(parents=True, exist_ok=True)
                jobs = [(fixtures, raw, binary[arm], q, r, expected, arm, CAL_SEED, cap)
                        for q, r, expected in sources]
                for row in pool.map(replay.dock, jobs):
                    calibration.append(dict(arm=arm, cap=cap, query=row["query"], receiver=row["receiver"], cpu_s=row["cpu_s"]))
            seeds = replay.SEEDS
        replay.reciprocal.write_csv(output / "calibration.csv", calibration)
        timings = defaultdict(lambda: defaultdict(float))
        for row in calibration:
            key = (row["query"], row["receiver"]) if args.per_complex_calibration else ("all", "all")
            timings[key][row["arm"], int(row["cap"])] += float(row["cpu_s"])
        caps, factors = {}, []
        for key, times in timings.items():
            slope = (times["native", HIGH] - times["native", LOW]) / (HIGH - LOW)
            if slope <= 0:
                raise ValueError("Positive search time increment required")
            cap = max(1000, round((LOW + (times["official", LOW] - times["native", LOW]) / slope) / 1000) * 1000)
            caps[key] = cap
            factors.append(dict(query=key[0], receiver=key[1], calibration_seed=CAL_SEED, cap=cap,
                               official_cpu_s_at_low=times["official", LOW], native_cpu_s_at_low=times["native",LOW],
                               native_cpu_s_at_high=times["native",HIGH]))
        serialized = factors if args.per_complex_calibration else {k:v for k,v in factors[0].items() if k not in ("query", "receiver")}
        (output / "frozen_cap.json").write_text(json.dumps(serialized, indent=2) + "\n")
        jobs = []
        raw = output / "raw"
        raw.mkdir(exist_ok=True)
        for i,(query, receiver, expected) in enumerate(sources):
            cap = caps[query["pdb"], receiver["pdb"]] if args.per_complex_calibration else caps["all", "all"]
            for j, seed in enumerate(seeds):
                arms = [("official", LOW), ("native", cap), ("geometry_analytic", cap)]
                for arm, budget in (arms if (i+j) % 2 == 0 else reversed(arms)):
                    jobs.append((fixtures, raw, binary[arm], query, receiver, expected, arm, seed, budget))
        rows = list(pool.map(replay.dock, jobs))
        replay.reciprocal.write_csv(output / "cells.csv", rows)
    result = {}
    for arm in binary:
        data = [r for r in rows if r["arm"] == arm]
        result[arm] = dict(n=len(data), top1=sum(r["first"]<=2 for r in data),
                           top5=sum(r["best5"]<=2 for r in data), oracle=sum(r["best"]<=2 for r in data),
                           models=sum(r["modes"] for r in data), cpu_s=round(sum(r["cpu_s"] for r in data),3))
    (output / "summary.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
