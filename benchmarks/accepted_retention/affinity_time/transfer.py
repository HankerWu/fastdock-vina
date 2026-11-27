"""Execute receiver-capped Gal3 arms using only archived calibration times."""

import argparse
import hashlib
import json
import os
import random
import sys
import time
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "affinity"))
from run import COUNT, cells, dock, read, write

__all__ = ["COUNT", "HERE", "METHODS", "SEED", "calibrate", "read", "write"]

SEED = 411119
METHODS = ("parent", "short_parent", "retention", "recurrence")


def calibrate():
    r = read(HERE.parent / "affinity/cost.csv")
    values = {
        (x["receiver"], x["method"]): Decimal(x["elapsed_s"])
        for x in r
        if x["panel"] == "Gal3"
    }
    rows = []
    for rec in ("1KJL", "3T1L"):
        parent = values[rec, "parent"]
        candidate = values[rec, "retention"]
        estimate = Decimal(12000) * parent / candidate
        cap = int(min(12000, estimate) // 500) * 500
        rows.append(
            {
                "receiver": rec,
                "parent_s": str(parent),
                "retention_s": str(candidate),
                "raw_projected_cap": str(estimate),
                "cap": cap,
                "saturated_at_12000": cap == 12000,
            }
        )
    return rows


def main(args):
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    (work / "raw").mkdir(exist_ok=args.full_cap)
    if args.full_cap and not (work / "runs.csv").exists():
        raise ValueError(
            "Run the primary four-arm protocol before the factorial completion."
        )
    os.sched_setaffinity(0, {args.cpu})
    start = time.perf_counter_ns()
    calibration = calibrate()
    overhead = time.perf_counter_ns() - start
    if not args.full_cap:
        write(work / "calibration.csv", calibration)
    caps = {r["receiver"]: r["cap"] for r in calibration}
    assert caps == {"1KJL": 11000, "3T1L": 11500}
    entries = cells(work, "Gal3")
    random.Random(SEED).shuffle(entries)
    binaries = {
        m: args.parent if m in ("parent", "short_parent") else getattr(args, m)
        for m in METHODS
    }
    methods = METHODS
    if args.full_cap:
        methods = ("full_retention", "full_recurrence")
        binaries = {m: getattr(args, m.removeprefix("full_")) for m in methods}
    rows = []
    for i, cell in enumerate(entries):
        j = i % len(methods)
        for m in methods[j:] + methods[:j]:
            cap = 12000 if args.full_cap or m == "parent" else caps[cell[1]]
            r = dock(cell, work, binaries[m], m, SEED, cap, True)
            if not r["joint"]:
                raise ValueError("unrecorded joint evaluation count")
            rows.append(r)
    write(work / ("full_cap_runs.csv" if args.full_cap else "runs.csv"), rows)
    (work / ("full_cap_setup.json" if args.full_cap else "setup.json")).write_text(
        json.dumps(
            {
                "cpu": args.cpu,
                "calibration_compute_ns": overhead,
                "calibration_sha256": hashlib.sha256(
                    (HERE.parent / "affinity/cost.csv").read_bytes()
                ).hexdigest(),
                "binaries": {
                    m: hashlib.sha256(Path(v).read_bytes()).hexdigest()
                    for m, v in binaries.items()
                },
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--cpu", type=int, default=7)
    p.add_argument(
        "--full-cap",
        action="store_true",
        help="Run the prespecified follow-up block at 12000.",
    )
    for arm in ("parent", "retention", "recurrence"):
        p.add_argument("--" + arm, required=True, type=Path)
    main(p.parse_args())
