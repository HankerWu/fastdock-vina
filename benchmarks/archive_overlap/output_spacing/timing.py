"""Counterbalanced, counter-free timings for all 94 validated four-trajectory conditions."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys

from compare import BENCH, output

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("inputs", "verified", "out", "old", "new"):
        p.add_argument(key, type=Path)
    p.add_argument("--cpu", type=int, required=True)
    a = p.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    a.out.mkdir(parents=True, exist_ok=True)
    verified = {
        (r["condition"], r["arm"]): r["sha256"]
        for r in csv.DictReader((a.verified / "cells.csv").open())
    }
    keys = sorted({k for k, _ in verified})
    rows = []
    for repeat in range(3):
        rng = random.Random(21531 + repeat)
        order = keys.copy()
        rng.shuffle(order)
        for i, key in enumerate(order):
            _, receiver, query, seed, _, cap = key.split("_")
            for arm in ("old", "new")[:: 1 if (i + repeat) % 2 == 0 else -1]:
                binary = getattr(a, arm)
                log = a.out / "last.log"
                usage = a.out / "last.json"
                pdb = a.out / "last.pdbqt"
                command = [
                    binary,
                    "--receptor",
                    a.inputs / (receiver + ".pdbqt"),
                    "--config",
                    a.inputs / (receiver + ".box"),
                    "--ligand",
                    a.inputs / (query + ".pdbqt"),
                    "--out",
                    pdb,
                    "--seed",
                    seed,
                    "--cpu",
                    "1",
                    "--exhaustiveness",
                    "4",
                    "--max_evals",
                    cap,
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--verbosity",
                    "0",
                ]
                subprocess.run(
                    [
                        sys.executable,
                        "-S",
                        BENCH / "rigid_flip/process_cost.py",
                        log,
                        usage,
                        *command,
                    ],
                    check=True,
                )
                if hashlib.sha256(pdb.read_bytes()).hexdigest() != verified[key, arm]:
                    raise ValueError("Output drift " + key)
                rows.append(
                    dict(
                        repeat=repeat + 1,
                        condition=key,
                        arm=arm,
                        **json.loads(usage.read_text()),
                    )
                )
        output(a.out / "repeats.csv", rows)
    (a.out / "builds.json").write_text(
        json.dumps(
            {
                k: hashlib.sha256(getattr(a, k).read_bytes()).hexdigest()
                for k in ["old", "new"]
            },
            indent=2,
        )
        + "\n"
    )
