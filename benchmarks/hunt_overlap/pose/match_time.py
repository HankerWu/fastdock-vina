"""Time-only NR1H2 calibration followed by a separate-seed, budget-capped replay."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
from statistics import median

from panel import CAPS, HOLDOUT, configure_arm, digest, graph_and_receptor_tools


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3])
    p.add_argument("--cpu", type=int, required=True)
    a = p.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    work = a.work.resolve()
    panel, _ = graph_and_receptor_tools(a.repo.resolve(), HOLDOUT)
    paired = panel.paired
    arms = ("parent", "fit", "invalid")
    meta = {q: json.loads((work / q / "preparation.json").read_text()) for q in HOLDOUT}
    timing = work / "matched"
    (timing / "runs").mkdir(parents=True, exist_ok=True)
    for q in HOLDOUT:
        if not (timing / q).exists():
            (timing / q).symlink_to(work / q, target_is_directory=True)
    with (work / "cells.csv").open() as handle:
        original = {r["key"]: r["pdbqt_sha256"] for r in csv.DictReader(handle)}
    trials, budgets = [], []

    def measure_pair(q, r):
        measurements = {}

        def run_block(seed, caps, phase):
            for repeat in range(3):
                # Latin rotation gives each arm each position once.
                for arm in arms[repeat:] + arms[:repeat]:
                    cap = caps[arm]
                    configure_arm(arm)
                    row, _ = paired.execute(
                        timing, meta, q, r, seed, cap, arm, a.binary, False
                    )
                    key = row["key"]
                    raw = digest(timing / "runs" / (key + ".pdbqt"))
                    if key in original and original[key] != raw:
                        raise ValueError("Timing replay changed predictions: " + key)
                    samples = measurements.setdefault((seed, arm, cap), [])
                    if samples and samples[0][1] != raw:
                        raise ValueError("Unstable replay: " + key)
                    samples.append((row["seconds"], raw))
                    trials.append(
                        {**row, "phase": phase, "repeat": repeat, "pdbqt_sha256": raw}
                    )

        def halfway(arm, cap):
            return median(v[0] for v in measurements[625231, arm, cap])

        for cap in CAPS:
            run_block(625231, {arm: cap for arm in arms}, "calibration")
        target = halfway("parent", CAPS[1])
        caps = {"parent": CAPS[1]}
        for arm in arms[1:]:
            lo, hi = [halfway(arm, b) for b in CAPS]
            if hi <= target:
                cap = CAPS[1]
            elif lo >= target or hi <= lo:
                cap = CAPS[0]
            else:
                cap = math.floor(
                    CAPS[0] + (CAPS[1] - CAPS[0]) * (target - lo) / (hi - lo)
                )
            caps[arm] = cap
            budgets.append(
                {
                    "query": q,
                    "receiver": r,
                    "arm": arm,
                    "low_seconds": lo,
                    "high_seconds": hi,
                    "parent_target_seconds": target,
                    "cap": cap,
                }
            )
        run_block(625237, caps, "assessment")

    for q in HOLDOUT:
        for r in HOLDOUT:
            measure_pair(q, r)
    paired.write(work / "matched_trials.csv", trials)
    paired.write(work / "matched_budgets.csv", budgets)


if __name__ == "__main__":
    main()
