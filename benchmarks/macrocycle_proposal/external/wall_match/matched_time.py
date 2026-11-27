"""Calibrate using time/evaluation counters only; then run the frozen validation."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
import re
import subprocess
import time
from pathlib import Path
from statistics import median

TARGETS = ("5T9W", "7MBO")
SEEDS = (458921, 740317, 191009)
CAPS = (2000, 3000, 4500, 6000, 7500, 9000, 10500, 12000)
CALIBRATION_SEED = 281291
COUNTERS = (
    "hunt_calls",
    "hunt_evals",
    "second_calls",
    "second_evals",
    "explicit_calls",
    "explicit_evals",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, rows):
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def run(args, target, arm, seed, cap, rep, slot, pair_cap):
    data = args.inputs / target
    # Only the known box centre is used. Crystal RMSDs are analyzed separately.
    center = json.loads((data / "preparation.json").read_text())["center"]
    box = [
        x
        for a, c in zip("xyz", center)
        for x in ("--center_" + a, str(c), "--size_" + a, "26")
    ]
    key = f"{target}_{seed}_{pair_cap}_{rep}_{arm}"
    raw = args.output / "raw"
    path = raw / (key + ".pdbqt")
    cmd = [
        getattr(args, arm),
        "--receptor",
        data / "receptor.pdbqt",
        "--ligand",
        data / "ligand.pdbqt",
        *box,
        "--seed",
        str(seed),
        "--cpu",
        "1",
        "--exhaustiveness",
        "1",
        "--num_modes",
        "5",
        "--energy_range",
        "100",
        "--max_evals",
        str(cap),
        "--verbosity",
        "0",
        "--search_counters",
        "--out",
        path,
    ]
    start = time.perf_counter_ns()
    log = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
    elapsed = time.perf_counter_ns() - start
    (raw / (key + ".log")).write_text(log)
    counters = {k: int(v) for k, v in re.findall(r"(\w+)=(\d+)", log) if k in COUNTERS}
    if len(counters) != len(COUNTERS):
        raise ValueError("Incomplete counters: " + key)
    # Hash only; calibration never parses energies or coordinates.
    return {
        "key": key,
        "pdb": target,
        "arm": arm,
        "seed": seed,
        "cap": cap,
        "pair_cap": pair_cap,
        "rep": rep,
        "slot": slot,
        "elapsed_ns": elapsed,
        **counters,
        "joint": sum(counters[k] for k in COUNTERS if k.endswith("evals")),
        "sha256": sha(path),
    }


def main(args):
    os.sched_setaffinity(0, {args.cpu})
    (args.output / "raw").mkdir(parents=True, exist_ok=False)
    binary_hashes = {a: sha(getattr(args, a)) for a in ("lock", "parent")}
    provenance = {
        "cpu": args.cpu,
        "binaries": binary_hashes,
        "protocol_sha256": sha(Path(__file__).with_name("PROTOCOL.md")),
        "inputs": {
            f"{t}/{f}": sha(args.inputs / t / f)
            for t in TARGETS
            for f in ("ligand.pdbqt", "receptor.pdbqt", "preparation.json")
        },
    }
    (args.output / "metadata.json").write_text(json.dumps(provenance, indent=2) + "\n")
    if args.mode == "calibrate":
        warm = [
            run(args, t, a, CALIBRATION_SEED, 12000, -1, i, 12000)
            for t in TARGETS
            for i, a in enumerate(("parent", "lock"))
        ]
        write(args.output / "warmup.csv", warm)
        cells = [(t, CALIBRATION_SEED, c) for t in TARGETS for c in CAPS]
    else:
        selection_path = args.selection / "selection.json"
        if (
            sha(selection_path)
            != (args.selection / "selection.sha256").read_text().strip()
        ):
            raise ValueError("Changed calibration decision")
        if json.loads((args.selection / "metadata.json").read_text()) != provenance:
            raise ValueError(
                "Validation input, hardware binding, binary or protocol differs"
            )
        sel = json.loads(selection_path.read_text())
        cells = [(t, s, sel[t]["cap"]) for t in TARGETS for s in SEEDS]
    rows = []
    rng = random.Random(281291)
    for rep in range(4):
        order = list(cells)
        rng.shuffle(order)
        for t, s, c in order:
            original_index = cells.index((t, s, c))
            first_lock = (rep + original_index) % 2 == 1
            arms = ("lock", "parent") if first_lock else ("parent", "lock")
            for i, a in enumerate(arms):
                rows.append(
                    run(args, t, a, s, 12000 if a == "parent" else c, rep, i, c)
                )
            write(args.output / "all_runs.csv", rows)
    for key, group in itertools.groupby(
        sorted(rows, key=lambda x: (x["pdb"], x["seed"], x["arm"], x["cap"])),
        key=lambda x: (x["pdb"], x["seed"], x["arm"], x["cap"]),
    ):
        group = list(group)
        if any(len({r[k] for r in group}) != 1 for k in (*COUNTERS, "sha256")):
            raise ValueError(f"Non-deterministic run {key}")
    if args.mode == "calibrate":
        selected = {}
        costs = []
        for t in TARGETS:
            table = []
            for c in CAPS:
                items = {
                    (r["rep"], r["arm"]): r
                    for r in rows
                    if r["pdb"] == t and r["pair_cap"] == c
                }
                ratio = [
                    items[i, "lock"]["elapsed_ns"] / items[i, "parent"]["elapsed_ns"]
                    for i in range(4)
                ]
                eval_ratio = items[0, "lock"]["joint"] / items[0, "parent"]["joint"]
                ok = (
                    median(ratio) <= 0.90
                    and sum(v <= 0.95 for v in ratio) >= 3
                    and eval_ratio <= 1
                )
                row = {
                    "pdb": t,
                    "cap": c,
                    **{f"wall_ratio_r{i}": v for i, v in enumerate(ratio)},
                    "median_wall_ratio": median(ratio),
                    "eval_ratio": eval_ratio,
                    "qualifies": ok,
                }
                costs.append(row)
                if ok:
                    table.append(c)
            selected[t] = {
                "cap": max(table) if table else min(CAPS),
                "qualified": bool(table),
            }
        write(args.output / "costs.csv", costs)
        decision = args.output / "selection.json"
        decision.write_text(json.dumps(selected, indent=2) + "\n")
        (args.output / "selection.sha256").write_text(sha(decision) + "\n")
    else:
        # Feed one representative geometry per identical repeated output to the standalone checker.
        result = []
        for r in rows:
            if r["rep"] == 0:
                n = (
                    (args.output / "raw" / (r["key"] + ".pdbqt"))
                    .read_bytes()
                    .count(b"ENDMDL")
                )
                result.append({**r, "modes": n})
        write(args.output / "runs.csv", result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["calibrate", "validate"])
    for k in ("inputs", "output", "parent", "lock"):
        p.add_argument(k, type=Path)
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    p.add_argument("--selection", type=Path)
    a = p.parse_args()
    if a.mode == "validate" and a.selection is None:
        p.error("Set --selection to the frozen calibration directory")
    main(a)
