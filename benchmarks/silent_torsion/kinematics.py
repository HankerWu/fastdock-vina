"""Check full-coordinate/trajectory equivalence and process cost in randomized paired blocks."""

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH / "rigid_flip"))
from holdout import args_for, paired


def read(src):
    with src.open() as stream:
        return list(csv.DictReader(stream))


def jobs(args):
    cubes = []
    for r in read(args.pose / "cells.csv"):
        if r["arm"] != "parent":
            continue
        q, rec = r["query"], r["receiver"]
        meta = {rec: json.loads((args.pose / rec / "preparation.json").read_text())}
        cmd = args_for(args.pose, meta, q, rec, r["seed"], r["cap"], Path("dummy"))
        idx = cmd.index("--out")
        del cmd[idx : idx + 2]
        cubes.append((r["key"], cmd, r))
    for r in read(args.affinity / "cells.csv"):
        if r["arm"] != "parent":
            continue
        rec, lig = r["receptor"], r["compound"]
        cmd = [
            "--config",
            str(args.inputs / (rec + ".box")),
            "--receptor",
            str(args.inputs / (rec + ".pdbqt")),
            "--ligand",
            str(args.inputs / (lig + ".pdbqt")),
            "--cpu",
            "1",
            "--exhaustiveness",
            "1",
            "--seed",
            "214039",
            "--max_evals",
            r["cap"],
            "--num_modes",
            "9",
            "--energy_range",
            "100",
        ]
        cubes.append((rec + "_" + lig, cmd, r))
    return cubes


def measure(args, cases, cpu):
    os.sched_setaffinity(0, {cpu})
    rows = []
    order = random.Random(917713 + cpu)
    for rep in range(1, 4):
        order.shuffle(cases)
        for key, cmd, ref in cases:
            arms = [("parent", args.parent), ("deferred", args.deferred)]
            order.shuffle(arms)
            for arm, executable in arms:
                stem = args.out / (key + "_" + arm)
                pose = stem.with_suffix(".pdbqt")
                log = stem.with_suffix(".log")
                cost = stem.with_suffix(".cost")
                subprocess.run(
                    [
                        sys.executable,
                        "-S",
                        str(BENCH / "rigid_flip/process_cost.py"),
                        str(log),
                        str(cost),
                        str(executable),
                        *cmd,
                        "--search_counters",
                        "--out",
                        str(pose),
                    ],
                    check=True,
                )
                digest = hashlib.sha256(pose.read_bytes()).hexdigest()
                counts = paired.counts(log.read_text())
                if digest != ref["sha256"] or any(
                    counts[k] != int(ref[k]) for k in counts
                ):
                    raise ValueError(
                        f"Coordinate or search parity failed in {rep} {key} {arm}"
                    )
                rows.append(
                    dict(
                        repeat=rep,
                        key=key,
                        arm=arm,
                        cpu=cpu,
                        sha256=digest,
                        **counts,
                        **json.loads(cost.read_text()),
                    )
                )
                pose.unlink()
                log.unlink()
                cost.unlink()
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("out", "pose", "affinity", "inputs", "parent", "deferred"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--cpus", type=int, nargs="+", required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)
    cases = jobs(args)
    random.Random(90311).shuffle(cases)
    with ProcessPoolExecutor(len(args.cpus)) as pool:
        futures = [
            pool.submit(measure, args, cases[i :: len(args.cpus)], cpu)
            for i, cpu in enumerate(args.cpus)
        ]
        paired.write(args.out / "timing.csv", [r for f in futures for r in f.result()])
