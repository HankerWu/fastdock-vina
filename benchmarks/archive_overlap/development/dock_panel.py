"""Dock historical PTP1B/ESR1/PARP1 pairs and verify full-precision observers."""

import argparse
import gzip
import hashlib
import itertools
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
sys.path.insert(0, str(BENCH / "rigid_flip"))
import holdout


def counts(text):
    trajectories = [
        dict(re.findall(r"(\w+)=(\d+)", line))
        for line in text.splitlines()
        if line.startswith("Search trajectory:")
    ]
    if [int(row["index"]) for row in trajectories] != [1, 2, 3, 4]:
        raise ValueError("Expected four trajectory counter records")
    result = holdout.paired.counts(text)
    for key in ("hunt_calls", "hunt_evals", "second_calls", "second_evals"):
        result[key] = sum(int(row[key]) for row in trajectories)
    return result


def run(args):
    work = args.inputs.resolve()
    dest = args.output.resolve()
    raw = dest / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    models, rows = [], []
    metadata = {
        q: json.loads((work / q / "preparation.json").read_text())
        for q in holdout.CASES
    }
    # The established receptor-frame evaluator finds each receptor next to the runs directory.
    for code in metadata:
        link = dest / code
        if not link.exists():
            link.symlink_to(work / code, target_is_directory=True)
    for q, r in itertools.product(sorted(metadata), repeat=2):
        if (
            metadata[r]["receptor_status"] != "ready"
            or holdout.CASES[q][0] != holdout.CASES[r][0]
        ):
            continue
        try:
            holdout.paired.base.alignment(metadata[q], metadata[r])
        except ValueError:
            continue
        for j, (seed, cap) in enumerate(itertools.product(holdout.SEEDS, holdout.CAPS)):
            for arm in ("parent", "fixed") if j % 2 == 0 else ("fixed", "parent"):
                key = f"{q}_{r}_{seed}_{cap}_{arm}"
                prefix = raw / key
                out = prefix.with_suffix(".pdbqt")
                cli = holdout.args_for(work, metadata, q, r, seed, cap, out)
                cli[cli.index("--exhaustiveness") + 1] = "4"
                cli.append("--search_counters")
                log = prefix.with_suffix(".log")
                cost = prefix.with_suffix(".cost")
                plain, observe = (
                    (args.parent, args.parent_obs)
                    if arm == "parent"
                    else (args.fixed, args.fixed_obs)
                )
                subprocess.run(
                    [
                        sys.executable,
                        "-S",
                        str(BENCH / "rigid_flip/process_cost.py"),
                        str(log),
                        str(cost),
                        str(plain.resolve()),
                        *cli,
                    ],
                    check=True,
                )
                trace = prefix.with_suffix(".trace")
                observed = prefix.with_suffix(".observed")
                cli[cli.index("--out") + 1] = str(observed)
                result = subprocess.run(
                    [str(observe.resolve()), *cli],
                    env={**os.environ, "FD_ARCHIVE_TRACE": str(trace)},
                    capture_output=True,
                    text=True,
                    check=True,
                )
                work_counts = counts(log.read_text())
                if out.read_bytes() != observed.read_bytes() or work_counts != counts(
                    result.stdout
                ):
                    raise ValueError(f"Observer changed poses or work in {key}")
                observed.unlink()
                with (
                    trace.open("rb") as src,
                    gzip.open(str(trace) + ".gz", "wb") as dst,
                ):
                    shutil.copyfileobj(src, dst)
                trace.unlink()
                row = dict(
                    key=key,
                    query=q,
                    receiver=r,
                    seed=seed,
                    cap=cap,
                    arm=arm,
                    sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
                    **work_counts,
                    **json.loads(cost.read_text()),
                )
                models.extend(
                    holdout.paired.inspect(raw, metadata[q], metadata[r], row)
                )
                rows.append(row)
                print(key, flush=True)
    holdout.paired.write(dest / "cells.csv", rows)
    holdout.paired.write(dest / "modes.csv", models)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "output", "parent", "fixed", "parent-obs", "fixed-obs"):
        p.add_argument("--" + name, type=Path, required=True)
    run(p.parse_args())
