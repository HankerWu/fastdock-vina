"""Per-compound CPU-affine wall-time calibration and independent-seed paired runs."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from assay_transfer import RECEIVERS, save

LO, HI = 3000, 12000


def dock(work, fixtures, binaries, rec, comp, seed, arm, cap):
    key = f"{rec}-{comp}-{seed}-{arm}-{cap}"
    out = work / "raw" / (key + ".pdbqt")
    cmd = [
        binaries[arm],
        "--receptor",
        fixtures / (rec + ".pdbqt"),
        "--ligand",
        fixtures / f"{comp}.pdbqt",
        "--config",
        fixtures / (rec + ".box"),
        "--out",
        out,
        "--seed",
        str(seed),
        "--cpu",
        "1",
        "--exhaustiveness",
        "1",
        "--max_evals",
        str(cap),
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--verbosity",
        "0",
        "--search_counters",
    ]
    start = time.perf_counter_ns()
    p = subprocess.run(cmd, capture_output=True, check=True, text=True)
    ns = time.perf_counter_ns() - start
    log = p.stdout + p.stderr
    out.with_suffix(".log").write_text(log)
    counters = {
        k: int(v)
        for k, v in re.findall(
            r"\b(hunt_evals|second_evals|explicit_evals|hunt_calls|second_calls|explicit_calls)=(\d+)",
            log,
        )
    }
    assert len(counters) == 6
    energies = re.findall(r"^REMARK VINA RESULT:\s+([-\d.]+)", out.read_text(), re.M)
    return dict(
        key=key,
        receptor=rec,
        compound=comp,
        seed=seed,
        arm=arm,
        cap=cap,
        elapsed_ns=ns,
        first_score=energies[0],
        modes=len(energies),
        **counters,
    )


def calibration(work, fixtures, binaries, cells):
    rows, budgets = [], []
    arms = list(itertools.product(("parent", "linear"), (LO, HI)))
    for index, (rec, comp) in enumerate(cells):
        pos = index % len(arms)
        times = {}
        for arm, cap in arms[pos:] + arms[:pos]:
            r = dock(work, fixtures, binaries, rec, comp, 75467, arm, cap)
            rows.append(r)
            times[arm, cap] = r["elapsed_ns"]
        low, high, target = (
            times[a, c] for a, c in (("linear", LO), ("linear", HI), ("parent", HI))
        )
        if high <= target:
            cap, reason = HI, "high_within_target"
        elif low >= target or low >= high:
            cap, reason = LO, "low_or_nonmonotonic"
        else:
            # Floor with integer division; all times are recorded as integer ns.
            cap, reason = LO + (HI - LO) * (target - low) // (high - low), "interpolate"
        budgets.append(
            dict(
                receptor=rec,
                compound=comp,
                parent_high_ns=target,
                linear_low_ns=low,
                linear_high_ns=high,
                calibrated_cap=cap,
                reason=reason,
            )
        )
        if (index + 1) % 13 == 0:
            print(f"calibration {index + 1}/{len(cells)}", flush=True)
    save(work / "calibration.csv", rows)
    save(work / "budgets.csv", budgets)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fixtures", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--parent", type=Path, required=True)
    p.add_argument("--linear", type=Path, required=True)
    p.add_argument("--cpu", type=int, default=7)
    args = p.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    if os.sched_getaffinity(0) != {args.cpu}:
        raise RuntimeError("CPU affinity not applied")
    work = args.work.resolve()
    (work / "raw").mkdir(parents=True, exist_ok=True)
    binaries = {arm: getattr(args, arm).resolve() for arm in ("parent", "linear")}
    cells = list(itertools.product(sum(RECEIVERS.values(), ()), range(5, 18)))
    random.Random(911251).shuffle(cells)
    start = time.perf_counter()
    calibration(work, args.fixtures.resolve(), binaries, cells)
    digest = hashlib.sha256((work / "budgets.csv").read_bytes()).hexdigest()
    print("frozen budgets SHA256", digest, flush=True)
    with (work / "budgets.csv").open() as stream:
        budgets = {
            (r["receptor"], int(r["compound"])): int(r["calibrated_cap"])
            for r in csv.DictReader(stream)
        }
    cells.sort()
    random.Random(187097).shuffle(cells)
    rows = []
    for i, (rec, comp) in enumerate(cells):
        for arm in ("parent", "linear") if i % 2 == 0 else ("linear", "parent"):
            cap = HI if arm == "parent" else budgets[rec, comp]
            rows.append(
                dock(
                    work, args.fixtures.resolve(), binaries, rec, comp, 98251, arm, cap
                )
            )
    save(work / "validation.csv", rows)
    (work / "run.json").write_text(
        json.dumps(
            dict(
                cpu=args.cpu,
                protocol_cpu_mask=sorted(os.sched_getaffinity(0)),
                budget_sha256=digest,
                elapsed_experiment_s=time.perf_counter() - start,
                binary_sha256={
                    arm: hashlib.sha256(v.read_bytes()).hexdigest()
                    for arm, v in binaries.items()
                },
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
