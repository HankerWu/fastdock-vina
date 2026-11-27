"""Run the frozen two-ligand human FKBP12 experiment; each arm uses one process per cell."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
import re
import subprocess
import sys
import time
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def execute(args):
    sys.path.insert(0, str(args.repo.resolve() / "benchmarks/macrocycle_labels"))
    from run_labels import REFINE, TRAJECTORY

    os.sched_setaffinity(0, {args.cpu})
    root = args.work.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    raw = out / "raw"
    raw.mkdir(exist_ok=True)
    arms = dict(
        upstream=args.upstream, parent=args.parent, ccd=args.ccd, reset=args.reset
    )
    if args.lock:
        arms = {"lock": args.lock}
    binaries = {**arms, "counts": args.counts}
    initial_hashes = {arm: digest(path) for arm, path in binaries.items()}
    cells = list(
        itertools.product(
            ("1FKJ", "1FKF"), ("1FKJ", "6VCU"), (740317, 191009), (3000, 12000)
        )
    )
    random.Random(740317).shuffle(cells)
    if args.lock:
        with (out / "runs.csv").open() as f:
            records = [r for r in csv.DictReader(f) if r["arm"] != "lock"]
    else:
        records = []
    for i, (rec, lig, seed, cap) in enumerate(cells):
        receptor, ligand = root / rec, root / lig
        box = receptor / "assay.box"
        center = json.loads((receptor / "preparation.json").read_text())["center"]
        box.write_text(
            "".join(f"center_{k} = {v}\nsize_{k} = 24\n" for k, v in zip("xyz", center))
        )
        keys = list(arms)
        shift = i % len(keys)
        for arm in keys[shift:] + keys[:shift]:
            key = f"{rec}_{lig}_{seed}_{cap}_{arm}"
            dest = raw / (key + ".pdbqt")
            command = [
                "--ligand",
                ligand / "ligand.pdbqt",
                "--receptor",
                receptor / "receptor.pdbqt",
                "--config",
                box,
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
            ]
            start = time.perf_counter_ns()
            log = subprocess.check_output(
                [
                    arms[arm],
                    *command,
                    "--out",
                    dest,
                    *([] if arm == "upstream" else ["--search_counters"]),
                ],
                stderr=subprocess.STDOUT,
                text=True,
            )
            ns = time.perf_counter_ns() - start
            (raw / (key + ".log")).write_text(log)
            if arm == "upstream":
                replay = raw / (key + ".replay")
                log = subprocess.check_output(
                    [args.counts, *command, "--out", replay],
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                if dest.read_bytes() != replay.read_bytes():
                    raise ValueError("Upstream counting build changed pool: " + key)
                replay.unlink()
                (raw / (key + ".counts")).write_text(log)
            t, e = TRAJECTORY.search(log), REFINE.search(log)
            scores = re.findall(
                r"^REMARK VINA RESULT:\s+([-+.\d]+)", dest.read_text(), re.M
            )
            records.append(
                dict(
                    key=key,
                    receiver=rec,
                    query=lig,
                    seed=seed,
                    cap=cap,
                    arm=arm,
                    top_score=scores[0],
                    modes=len(scores),
                    wall_ns=ns,
                    sha256=digest(dest),
                    hunt_calls=t[1],
                    hunt_evals=t[2],
                    second_calls=t[3],
                    second_evals=t[4],
                    explicit_calls=e[1],
                    explicit_evals=e[2],
                    joint=int(t[2]) + int(t[4]) + int(e[2]),
                )
            )
    if initial_hashes != {arm: digest(path) for arm, path in binaries.items()}:
        raise RuntimeError(
            "A binary changed during the experiment. Repeat with stable builds."
        )
    with (out / "runs.csv").open("w") as f:
        w = csv.DictWriter(f, records[0], lineterminator="\n")
        w.writeheader()
        w.writerows(records)
    binary_file = out / "binaries.json"
    hashes = json.loads(binary_file.read_text()) if args.lock else {}
    hashes.update(
        {arm: digest(p) for arm, p in {**arms, "counts": args.counts}.items()}
    )
    binary_file.write_text(json.dumps(hashes, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("repo", "work", "out", "upstream", "counts", "parent", "ccd", "reset"):
        p.add_argument(k, type=Path)
    p.add_argument(
        "--lock",
        type=Path,
        help="Append the subsequent fixed-ring arm to an existing experiment.",
    )
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    execute(p.parse_args())
