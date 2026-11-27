"""Separate within-trajectory maintenance from later aggregation with matched inputs and costs."""

import argparse
import csv
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path

from dock_panel import BENCH, counts, holdout


def main(args):
    directory = args.output.resolve()
    raw = directory / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    work = args.inputs.resolve()
    meta = {
        q: json.loads((work / q / "preparation.json").read_text())
        for q in holdout.CASES
    }
    for q in meta:
        if not (directory / q).exists():
            (directory / q).symlink_to(work / q, target_is_directory=True)
    cases = [r for r in csv.DictReader(args.cells.open()) if r["arm"] == "parent"]
    rng = random.Random(97281)
    rows, poses, costs = [], [], []
    for rep in range(args.repeats):
        rng.shuffle(cases)
        for case in cases:
            q, r = case["query"], case["receiver"]
            arms = ["parent", "save", "merge", "fixed"]
            rng.shuffle(arms)
            for arm in arms:
                key = case["key"].rsplit("_", 1)[0] + "_" + arm
                path = raw / (key + ".pdbqt")
                log = raw / (key + ".log")
                cost = raw / (key + ".cost")
                cli = holdout.args_for(
                    work, meta, q, r, case["seed"], case["cap"], path
                )
                cli[cli.index("--exhaustiveness") + 1] = "4"
                cli.append("--search_counters")
                previous = path.read_bytes() if rep else None
                subprocess.run(
                    [
                        sys.executable,
                        "-S",
                        str(BENCH / "rigid_flip/process_cost.py"),
                        str(log),
                        str(cost),
                        str(getattr(args, arm).resolve()),
                        *cli,
                    ],
                    check=True,
                )
                if rep and previous != path.read_bytes():
                    raise ValueError("Non-reproducible poses in " + key)
                costs.append(
                    dict(key=key, repeat=rep + 1, **json.loads(cost.read_text()))
                )
                if rep:
                    continue
                row = dict(
                    key=key,
                    arm=arm,
                    query=q,
                    receiver=r,
                    seed=case["seed"],
                    cap=case["cap"],
                    sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    **counts(log.read_text()),
                )
                poses.extend(holdout.paired.inspect(raw, meta[q], meta[r], row))
                rows.append(row)
    holdout.paired.write(directory / "cells.csv", rows)
    holdout.paired.write(directory / "modes.csv", poses)
    holdout.paired.write(directory / "timings.csv", costs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "cells", "output", "parent", "save", "merge", "fixed"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--repeats", type=int, default=3)
    main(p.parse_args())
