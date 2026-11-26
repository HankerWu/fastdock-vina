"""Calibrate an end-to-end wall-time-matched native search on a pinned logical CPU."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import resource
import subprocess
import time
from pathlib import Path
from statistics import median

CALIBRATION_LIGANDS = (
    "5",
    "1c",
    "4",
)  # Smallest, middle, largest (26/27/29 heavy atoms).
RECEIVERS = ("3D49", "2ZFF", "2ZC9")
POSE_LIGANDS = ("5", "1a", "1b", "3a", "4")
CONFIRMATION_SEEDS = (88091, 302123)
CALIBRATION_SEED = 101921
BASE_CAP = 50000
PROBE_CAP = 75000


def write(path, rows):
    with path.open("w") as h:
        w = csv.DictWriter(h, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def run(binary, inputs, work, label, rec, lig, seed, kind, cap, phase, slot):
    input_dir = inputs / "prepared"
    path = work / (label + ".pdbqt")
    cmd = [
        str(binary),
        "--receptor",
        str(input_dir / (rec + ".pdbqt")),
        "--ligand",
        str(input_dir / (lig + ".pdbqt")),
        "--config",
        str(input_dir / (rec + ".box")),
        "--cpu",
        "1",
        "--pose_selection",
        "late",
        "--exhaustiveness",
        "4",
        "--max_evals",
        str(cap),
        "--num_modes",
        "12",
        "--refine_candidates",
        "16",
        "--energy_range",
        "8",
        "--search_counters",
        "--seed",
        str(seed),
        "--out",
        str(path),
    ]
    if kind == "soft":
        cmd.append("--soft_receptor_search")
    # Parent and all children inherit the affinity mask selected in main().
    moment = time.perf_counter_ns()
    prev = resource.getrusage(resource.RUSAGE_CHILDREN)
    with path.with_suffix(".log").open("w") as h:
        subprocess.run(cmd, stdout=h, stderr=subprocess.STDOUT, check=True)
    elapsed = (time.perf_counter_ns() - moment) / 1e9
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "phase": phase,
        "slot": slot,
        "receipt": label,
        "receiver": rec,
        "ligand": lig,
        "seed": seed,
        "arm": kind,
        "cap": cap,
        "wall_seconds": elapsed,
        "child_cpu_seconds": after.ru_utime
        - prev.ru_utime
        + after.ru_stime
        - prev.ru_stime,
        "output_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def calibrate(rows, metric):
    grouped = {}
    for row in rows:
        grouped.setdefault((row["receiver"], row["ligand"], row["slot"]), {})[
            (row["arm"], row["cap"])
        ] = row[metric]
    delta_soft = [v["soft", BASE_CAP] - v["hard", BASE_CAP] for v in grouped.values()]
    delta_probe = [v["hard", PROBE_CAP] - v["hard", BASE_CAP] for v in grouped.values()]
    if sum(delta_probe) <= 0:
        raise RuntimeError("Native budget probe did not increase mean cost")
    increase = (PROBE_CAP - BASE_CAP) * sum(delta_soft) / sum(delta_probe)
    return {
        "extra_cap": increase,
        "estimated_cap": max(BASE_CAP, round((BASE_CAP + increase) / 1000) * 1000),
        "hard": sum(v["hard", BASE_CAP] for v in grouped.values()),
        "soft": sum(v["soft", BASE_CAP] for v in grouped.values()),
        "probe": sum(v["hard", PROBE_CAP] for v in grouped.values()),
        "median_paired_overhead": median(
            v["soft", BASE_CAP] / v["hard", BASE_CAP] - 1 for v in grouped.values()
        ),
    }


def main(args):
    w = args.work.resolve()
    w.mkdir(exist_ok=True, parents=True)
    bin = args.binary.resolve()
    inputs = args.inputs.resolve()
    allowed = os.sched_getaffinity(0)
    core = args.core if args.core is not None else max(allowed)
    os.sched_setaffinity(0, {core})
    assert os.sched_getaffinity(0) == {core}
    common = {"binary": bin, "inputs": inputs, "work": w}
    # Three discarded warm-ups cover both grid implementations and all three receptor inputs.
    for k, r in enumerate(RECEIVERS):
        run(
            **common,
            label=f"warmup_{r}",
            rec=r,
            lig=CALIBRATION_LIGANDS[k],
            seed=101917,
            kind="soft" if k == 1 else "hard",
            cap=BASE_CAP,
            phase="warmup",
            slot=k,
        )
    rows = []
    specs = [("hard", BASE_CAP), ("soft", BASE_CAP), ("hard", PROBE_CAP)]
    cells = list(itertools.product(RECEIVERS, CALIBRATION_LIGANDS))
    for block in (0, 1):
        for j, (r, lig) in enumerate(cells if block == 0 else reversed(cells)):
            i = (j + block) % 3
            ordered = specs[i:] + specs[:i]
            if block:
                ordered = list(reversed(ordered))
            for arm, cap in ordered:
                rows.append(
                    run(
                        **common,
                        label=f"cal_{r}_{lig}_{block}_{arm}_{cap}",
                        rec=r,
                        lig=lig,
                        seed=CALIBRATION_SEED,
                        kind=arm,
                        cap=cap,
                        phase="calibration",
                        slot=block,
                    )
                )
    write(w / "calibration.csv", rows)
    wall = calibrate(rows, "wall_seconds")
    cpu = calibrate(rows, "child_cpu_seconds")
    cfg = {
        "logical_cpu": core,
        "original_affinity": sorted(allowed),
        "binary_hash": hashlib.sha256(bin.read_bytes()).hexdigest(),
        "wall": wall,
        "cpu": cpu,
    }
    (w / "calibrated.json").write_text(json.dumps(cfg, indent=2) + "\n")
    print(json.dumps(cfg), flush=True)
    cap = wall["estimated_cap"]
    rr = []
    # A matched standard-cap ablation provides a same-core time check on the actual evaluation panel.
    arms = [
        ("hard", BASE_CAP, "base"),
        ("soft", BASE_CAP, "soft"),
        ("hard", cap, "hard_equal"),
    ]
    for j, (r, lig, seed) in enumerate(
        itertools.product(RECEIVERS, POSE_LIGANDS, CONFIRMATION_SEEDS)
    ):
        i = j % 3
        order = arms[i:] + arms[:i]
        if (j // 3) % 2:
            order = list(reversed(order))
        for method, budget, label in order:
            rr.append(
                run(
                    **common,
                    label=f"test_{r}_{lig}_{seed}_{label}",
                    rec=r,
                    lig=lig,
                    seed=seed,
                    kind=method,
                    cap=budget,
                    phase=label,
                    slot=j,
                )
            )
    write(w / "validation.csv", rr)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for n in ("binary", "inputs", "work"):
        p.add_argument(n, type=Path)
    p.add_argument("--core", type=int)
    main(p.parse_args())
