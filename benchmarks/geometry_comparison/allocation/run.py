"""Fixed total nominal search effort across breadth and depth on the human panel."""

import argparse
import itertools
import json
from pathlib import Path
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pose_run import align_with_numbering, load_meta, measure, write_csv

# Locked before running. New seeds, identical to both allocations.
SEEDS = (76231, 76243, 76249)
ARMS = (("depth", 3, 48000), ("breadth", 12, 12000))


def run(spec):
    task, trajectories = spec
    row, modes = measure(task, trajectories)
    row["trajectories"] = trajectories
    return row, modes


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fixtures", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    work = a.fixtures.resolve()
    manifest = json.loads((work / "manifest.json").read_text())
    metas = [
        load_meta(work, row["pdb"]) for row in manifest if row["status"] == "ready"
    ]
    # Reciprocal directions require a verified local backbone alignment.
    specs = []
    for q, r in itertools.product(metas, repeat=2):
        if q["target"] != r["target"]:
            continue
        rot, tr, rms, n, offset = np.eye(3), np.zeros(3), 0.0, 0, 0
        if q["pdb"] != r["pdb"]:
            rot, tr, rms, n, offset = align_with_numbering(q, r)
        ref = np.asarray(q["mapped"]) @ rot + tr
        for i, seed in enumerate(SEEDS):
            for arm, trajectories, cap in ARMS[:: 1 if i % 2 == 0 else -1]:
                task = (
                    work,
                    a.work.resolve(),
                    a.binary.resolve(),
                    q,
                    r,
                    ref,
                    rms,
                    n,
                    offset,
                    seed,
                    arm,
                    cap,
                )
                specs.append((task, trajectories))
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        results = list(pool.map(run, specs))
    write_csv(a.work / "cells.csv", [row for row, _ in results])
    write_csv(a.work / "modes.csv", [m for _, modes in results for m in modes])


if __name__ == "__main__":
    main()
