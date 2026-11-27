"""Replay the unchanged 32 Pro/CCD jobs for time only, reversing arm order on alternate rounds."""

import argparse
import csv
import hashlib
import itertools
import os
import random
import subprocess
import sys
import tempfile
import time
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
for k in ("repo", "work", "out", "parent", "ccd"):
    p.add_argument(k, type=Path)
p.add_argument("--lock", type=Path, help="Time Pro against fixed ring separately.")
p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
a = p.parse_args()
os.sched_setaffinity(0, {a.cpu})
sys.path.insert(0, str(a.repo.resolve() / "benchmarks/macrocycle_labels"))
from run_labels import REFINE, TRAJECTORY

with (a.out / "runs.csv").open() as f:
    src = {r["key"]: r for r in csv.DictReader(f)}
rows = []
with tempfile.TemporaryDirectory() as td:
    tmp = Path(td) / "pose.pdbqt"
    cells = list(
        itertools.product(
            ("1FKF", "1FKJ"), ("1FKJ", "6VCU"), (740317, 191009), (3000, 12000)
        )
    )
    random.Random(191009).shuffle(cells)
    for rep in range(4):
        for i, (rec, q, seed, cap) in enumerate(cells):
            arm_name = "lock" if a.lock else "ccd"
            order = ("parent", arm_name) if (i + rep) % 2 == 0 else (arm_name, "parent")
            for arm in order:
                key = f"{rec}_{q}_{seed}_{cap}_{arm}"
                old = src[key]
                start = time.perf_counter_ns()
                log = subprocess.check_output(
                    [
                        getattr(a, arm),
                        "--receptor",
                        a.work / rec / "receptor.pdbqt",
                        "--ligand",
                        a.work / q / "ligand.pdbqt",
                        "--config",
                        a.work / rec / "assay.box",
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
                        "--out",
                        tmp,
                        "--search_counters",
                    ],
                    text=True,
                    stderr=subprocess.STDOUT,
                )
                elapsed = time.perf_counter_ns() - start
                t, e = TRAJECTORY.search(log), REFINE.search(log)
                expected = [
                    old[v]
                    for v in (
                        "hunt_calls",
                        "hunt_evals",
                        "second_calls",
                        "second_evals",
                        "explicit_calls",
                        "explicit_evals",
                    )
                ]
                sha = hashlib.sha256(tmp.read_bytes()).hexdigest()
                if expected != [*t.groups(), *e.groups()] or sha != old["sha256"]:
                    raise RuntimeError("Nonreproducing output or counters: " + key)
                rows.append(
                    dict(
                        round=rep + 1,
                        key=key,
                        arm=arm,
                        query=q,
                        receiver=rec,
                        seed=seed,
                        cap=cap,
                        sha256=sha,
                        joint=old["joint"],
                        wall_ns=elapsed,
                    )
                )
with (a.out / ("lock_timing.csv" if a.lock else "repeated_timing.csv")).open("w") as f:
    writer = csv.DictWriter(f, rows[0], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
