"""Repeat four-trajectory timings without changing candidates or scoring the repeats."""

import argparse
import os
import random
import statistics
from pathlib import Path

from replay import cpu, digest
from run import cells, dock, read, write


def measure(args):
    os.sched_setaffinity(0, {args.cpu})
    source = read(args.work / "runs.csv")
    expected = {(row["receptor"], row["compound"], row["arm"]): row for row in source}
    entries = cells(args.work, "Gal3")
    rows = []
    for repeat in (1, 2):
        work = args.scratch / str(repeat)
        (work / "raw").mkdir(parents=True, exist_ok=True)
        random.Random(3100 + repeat).shuffle(entries)
        for index, cell in enumerate(entries):
            arms = (
                ("parent", "candidate")
                if (index + repeat) % 2
                else ("candidate", "parent")
            )
            for arm in arms:
                ref = expected[cell[1], cell[2], arm]
                before = cpu()
                row = dock(
                    cell,
                    work,
                    getattr(args, arm),
                    arm,
                    391781,
                    3000,
                    True,
                    trajectories=4,
                )
                usage = cpu() - before
                raw = (work / "raw" / (row["key"] + ".pdbqt")).read_bytes()
                assert digest(raw) == ref["pool_sha256"] and int(row["joint"]) == int(
                    ref["joint"]
                )
                rows.append(
                    {
                        "key": row["key"],
                        "repeat": repeat,
                        "arm": arm,
                        "elapsed_ns": row["elapsed_ns"],
                        "cpu_s": usage,
                    }
                )
    write(args.work / "timing.csv", rows)
    medians = []
    for arm in ("parent", "candidate"):
        selected = [row for row in source if row["arm"] == arm]
        wall_s, cpu_s = 0.0, 0.0
        for initial in selected:
            obs = [initial] + [r for r in rows if r["key"] == initial["key"]]
            assert len(obs) == 3
            wall_s += statistics.median(float(r["elapsed_ns"]) for r in obs) / 1e9
            cpu_s += statistics.median(float(r["cpu_s"]) for r in obs)
        medians.append(
            {
                "arm": arm,
                "n": len(selected),
                "repeats": 3,
                "sum_median_wall_s": wall_s,
                "sum_median_cpu_s": cpu_s,
            }
        )
    write(args.work / "timing_medians.csv", medians)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "scratch", "parent", "candidate"):
        p.add_argument(name, type=Path)
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    measure(p.parse_args())
