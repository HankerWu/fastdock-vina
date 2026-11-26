"""Repeat six dockers on identical cognate or reciprocal inputs and independently seeded runs."""
import argparse
import csv
import hashlib
import itertools
import json
import os
import re
import statistics
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "cubic_grid" / "confirmation"))
import run as reciprocal

SEEDS = (85231, 85237)
ARMS = ("official", "native", "late", "armijo", "geometry_fd", "geometry_analytic")


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def experiments(protocol, fixtures):
    if protocol == "cognate":
        # Explicit manifest avoids incorporating stale files from failed preparation.
        for entry in json.loads((fixtures / "manifest.json").read_text()):
            if entry["status"] == "ready":
                meta = json.loads((fixtures / entry["pdb"] / "metadata.json").read_text())
                yield meta, meta, np.asarray(meta["mapped"])
    else:
        data = [json.loads(path.read_text()) for path in sorted(fixtures.glob("*/preparation.json"))]
        for target in reciprocal.CORE_TARGETS:
            pair = [row for row in data if row["target"] == target and row["status"] == "ready"]
            if len(pair) != 2:
                raise ValueError(f"Expected two prepared receptors for {target}")
            for query, receiver in itertools.permutations(pair):
                yield query, receiver, reciprocal.reference(query, receiver)


def dock(job):
    fixtures, raw, binary, query, receiver, expected, arm, seed, cap = job
    key = f"{query['pdb']}-{receiver['pdb']}-{seed}-{arm}"
    pdbqt, log = raw / (key + ".pdbqt"), raw / (key + ".log")
    ligand = fixtures / query["pdb"] / "ligand.pdbqt"
    receptor = fixtures / receiver["pdb"] / "receptor.pdbqt"
    cmd = [str(binary), "--ligand", str(ligand), "--receptor", str(receptor), "--out", str(pdbqt),
           "--cpu", "1", "--seed", str(seed), "--max_evals", str(cap), "--exhaustiveness", "3",
           "--num_modes", "9", "--energy_range", "100", "--spacing", "0.375"]
    for axis, center in zip("xyz", receiver["center"]):
        cmd.extend(["--size_" + axis, "24", "--center_" + axis, str(center)])
    if arm != "official":
        cmd.append("--search_counters")
    if arm not in ("official", "native"):
        cmd.extend(["--pose_selection", "late"])
    if arm in ("armijo", "geometry_fd", "geometry_analytic"):
        cmd.extend(["--refinement_optimizer", "armijo" if arm == "armijo" else "geometry"])
    clock0, cpu0 = time.monotonic(), os.times()
    with log.open("w") as stream:
        subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=180)
    cpu1 = os.times()
    wall = time.monotonic() - clock0
    cpu = cpu1.children_user + cpu1.children_system - cpu0.children_user - cpu0.children_system
    content = log.read_text()
    rmsds = reciprocal.measure(pdbqt, query["serial"], expected)
    scores = [float(l.split()[3]) for l in pdbqt.read_text().splitlines() if l.startswith("REMARK VINA RESULT:")]
    if len(scores) != len(rmsds):
        raise ValueError(f"Pose/score count mismatch: {key}")
    hunt = re.findall(r"Search trajectory: .*?hunt_evals=(\d+).*?second_evals=(\d+)", content)
    counts = re.search(r"Search counters: .*?explicit_inputs=(\d+) explicit_calls=(\d+) explicit_evals=(\d+) late_prescore=(\d+)", content)
    if arm != "official" and (len(hunt) != 3 or not counts):
        raise ValueError(f"Incomplete accounting: {key}")
    return dict(target=query["target"], query=query["pdb"], receiver=receiver["pdb"], seed=seed, arm=arm,
                first=rmsds[0], best5=min(rmsds[:5]), best=min(rmsds), native_score=scores[0], modes=len(scores),
                hunt_evals=sum(int(h) for h,s in hunt) if hunt else None,
                second_evals=sum(int(s) for h,s in hunt) if hunt else None,
                explicit_inputs=int(counts[1]) if counts else None, explicit_calls=int(counts[2]) if counts else None,
                explicit_evals=int(counts[3]) if counts else None, prescore=int(counts[4]) if counts else None,
                cpu_s=round(cpu, 3), wall_s=round(wall,4), pool_sha256=checksum(pdbqt))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("protocol", choices=("cognate", "reciprocal"))
    for arg in ("fixtures", "output", "official", "branch", "analytic"):
        p.add_argument(arg, type=Path)
    p.add_argument("--jobs", type=int, default=1)
    args = p.parse_args()
    fixtures, output = args.fixtures.resolve(), args.output.resolve()
    binaries = {a: (args.official if a == "official" else args.analytic if a == "geometry_analytic" else args.branch).resolve() for a in ARMS}
    raw = output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    jobs = []
    for i, (query, receiver, expected) in enumerate(experiments(args.protocol, fixtures)):
        for j, seed in enumerate(SEEDS):
            for arm in (ARMS if (i+j) % 2 == 0 else reversed(ARMS)):
                jobs.append((fixtures, raw, binaries[arm], query, receiver, expected, arm, seed, 12000))
    if not jobs:
        raise ValueError("No eligible cells")
    with ProcessPoolExecutor(args.jobs) as pool:
        rows = list(pool.map(dock, jobs))
    reciprocal.write_csv(output / "cells.csv", rows)
    arms = {}
    for arm in ARMS:
        selected = [r for r in rows if r["arm"] == arm]
        arms[arm] = dict(n=len(selected), model_count=sum(r["modes"] for r in selected),
                         top1=sum(r["first"] <= 2 for r in selected), top5=sum(r["best5"] <= 2 for r in selected),
                         oracle=sum(r["best"] <= 2 for r in selected), median_top1=statistics.median(r["first"] for r in selected),
                         binary_sha256=checksum(binaries[arm]), child_cpu_s=round(sum(r["cpu_s"] for r in selected),3))
    (output / "summary.json").write_text(json.dumps(arms, indent=2) + "\n")


if __name__ == "__main__":
    main()
