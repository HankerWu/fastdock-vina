"""Replay frozen comparisons, aggregate all four trajectories and audit CLI outputs."""

import argparse
import hashlib
import itertools
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

from prepare_inputs import BENCH, SYSTEMS

sys.path.insert(0, str(BENCH / "pair_gradient/pose"))
import paired

ARMS = ("upstream", "parent", "candidate")
SEEDS = (684173, 961751)
CAP = 12000
TRAJECTORIES = 4


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cli(binary, common, out, log, counters=False):
    # Use a lightweight supervisor to avoid inheriting the analysis process RSS.
    argv = [str(binary), *common, "--out", str(out)]
    if counters:
        argv.append("--search_counters")
    stats = log.with_suffix(log.suffix + ".usage")
    subprocess.run(
        [
            sys.executable,
            "-S",
            str(BENCH / "rigid_flip/process_cost.py"),
            str(log),
            str(stats),
            *argv,
        ],
        check=True,
    )
    usage = json.loads(stats.read_text())
    return usage["wall_s"], usage["cpu_s"], usage["peak_rss_KiB"]


def sum_counts(log):
    lines = [line for line in log.splitlines() if "Search trajectory:" in line]
    if len(lines) != TRAJECTORIES:
        raise ValueError(f"Expected {TRAJECTORIES} trajectories, observed {len(lines)}")
    fields = [dict(re.findall(r"(\w+)=(\d+)", line)) for line in lines]
    out = {
        k: sum(int(row[k]) for row in fields)
        for k in ("hunt_calls", "hunt_evals", "second_calls", "second_evals")
    }
    overall = [line for line in log.splitlines() if "Search counters:" in line]
    if len(overall) != 1:
        raise ValueError("Missing explicit-stage counters")
    explicit = dict(re.findall(r"(\w+)=(\d+)", overall[0]))
    out.update(
        {
            k: int(explicit[k])
            for k in ("explicit_inputs", "explicit_calls", "explicit_evals")
        }
    )
    return out


def pool_spacing(path, serial):
    xyz = []
    for block in path.read_text().split("ENDMDL")[:-1]:
        atoms = paired.atoms(block.splitlines())
        xyz.append(np.array([atoms[key][0] for key in serial]))
    distances = [
        np.sqrt(np.square(a - b).sum(axis=1).mean())
        for a, b in itertools.combinations(xyz, 2)
    ]
    return (
        min(distances) if distances else None,
        sum(d < 1.0 - 0.002 for d in distances),
    )


def run(args, systems=SYSTEMS, seeds=SEEDS, arms=ARMS, queries=None):
    work = args.work.resolve()
    os.sched_setaffinity(0, {args.cpu})
    bins = {arm: getattr(args, arm).resolve() for arm in (*arms, "counter")}
    if "hold" in arms:
        bins["hold_counter"] = args.hold_counter.resolve()
    (work / "binaries.json").write_text(
        json.dumps({k: sha(v) for k, v in bins.items()}, indent=2) + "\n"
    )
    runs = work / "runs"
    runs.mkdir(exist_ok=True)
    meta = {
        code: json.loads((work / code / "preparation.json").read_text())
        for code in systems
        if (work / code / "preparation.json").is_file()
    }
    cells = [
        (q, r, s)
        for q, r in itertools.product(meta, repeat=2)
        if (queries is None or q in queries)
        and meta[q]["target"] == meta[r]["target"]
        and meta[r]["receptor_status"] == "ready"
        for s in seeds
    ]
    random.Random(seeds[0]).shuffle(cells)
    if args.reverse:
        cells.reverse()
    if args.rescore:
        old = json.loads((work / "timings.json").read_text())
    else:
        old = {}
    rows, modes = [], []
    for i, (q, r, seed) in enumerate(cells):
        common = [
            "--receptor",
            str(work / r / "receptor.pdbqt"),
            "--ligand",
            str(work / q / "ligand.pdbqt"),
            "--cpu",
            "1",
            "--exhaustiveness",
            str(TRAJECTORIES),
            "--max_evals",
            str(CAP),
            "--num_modes",
            "9",
            "--min_rmsd",
            "1.0",
            "--energy_range",
            "100",
            "--seed",
            str(seed),
        ]
        for axis, value in zip("xyz", meta[r]["center"]):
            common += ["--center_" + axis, str(value), "--size_" + axis, "22.5"]
        for arm in arms[i % len(arms) :] + arms[: i % len(arms)]:
            key = f"{q}_{r}_{seed}_{arm}"
            output = runs / (key + ".pdbqt")
            counts = runs / (key + ".counts")
            if not args.rescore:
                old[key] = cli(bins[arm], common, output, runs / (key + ".log"))
                replay = runs / (key + ".replay")
                cli(
                    bins[{"upstream": "counter", "hold": "hold_counter"}.get(arm, arm)],
                    common,
                    replay,
                    counts,
                    counters=arm != "upstream",
                )
                if sha(output) != sha(replay):
                    raise ValueError("Counted output not identical: " + key)
                replay.unlink()
                print(key, old[key], flush=True)
                (work / "timings.json").write_text(json.dumps(old, indent=2) + "\n")
            wall, cpu, rss = old[key]
            row = dict(
                key=key,
                target=meta[q]["target"],
                query=q,
                receiver=r,
                seed=seed,
                arm=arm,
                wall_s=wall,
                cpu_s=cpu,
                maxrss_KiB=rss,
                output_sha256=sha(output),
                **sum_counts(counts.read_text()),
            )
            if "hold" in arms:
                skips = [
                    int(n)
                    for n in re.findall(r"Skipped refills: (\d+)", counts.read_text())
                ]
                if arm == "hold" and len(skips) != TRAJECTORIES:
                    raise ValueError(f"Missing hold observer: {key}")
                row["skipped_refills"] = sum(skips) if skips else ""
                row["trajectories_skipping_refills"] = (
                    sum(n > 0 for n in skips) if skips else ""
                )
            poses = paired.inspect(runs, meta[q], meta[r], row)
            row["pool_min_indexed_rmsd_A"], row["pool_pairs_below_0.998_A"] = (
                pool_spacing(output, meta[q]["serial"])
            )
            rows.append(row)
            modes.extend(poses)
        paired.write(work / "cells.csv", rows)
        paired.write(work / "modes.csv", modes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    for arm in (*ARMS, "counter"):
        parser.add_argument("--" + arm, type=Path, required=True)
    parser.add_argument("--cpu", type=int, required=True)
    parser.add_argument("--rescore", action="store_true")
    parser.add_argument(
        "--reverse",
        action="store_true",
        help="Reverse the frozen cell order for timing repeats",
    )
    run(parser.parse_args())
