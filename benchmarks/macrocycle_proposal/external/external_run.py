"""Run the macrocycle proposal on every eligible external PDB accession."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
import re
import subprocess
import time
from pathlib import Path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, records):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=records[0], lineterminator="\n")
        w.writeheader()
        w.writerows(records)


def run(args):
    os.sched_setaffinity(0, {args.cpu})
    arms = {
        name: getattr(args, name).resolve()
        for name in ("parent", "ccd", "reset", "upstream")
    }
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    raw = output / "raw"
    raw.mkdir(exist_ok=True)
    if args.ids:
        ids = args.ids
    else:
        ids = sorted(
            p.parent.name
            for p in args.work.glob("*/screening.json")
            if json.loads(p.read_text())["status"] == "ready"
        )
    if not ids:
        raise ValueError("No eligible external structures")
    seeds = [458921, 740317, 191009]
    if args.mode == "persistent":
        if args.lock is None:
            raise ValueError("The persistent control requires --lock")
        arms = {"parent": arms["parent"], "lock": args.lock.resolve()}
        cells = list(itertools.product(ids, seeds, [2000, 3000, 9000, 12000]))
    elif args.mode == "paired":
        cells = list(itertools.product(ids, seeds, [3000, 12000]))
    elif args.mode == "lower":
        arms = {k: v for k, v in arms.items() if k in ("reset", "ccd")}
        cells = list(itertools.product(ids, seeds, [2000, 9000]))
    else:
        cells = list(itertools.product(ids, seeds, [3000, 12000]))
    random.Random(20261127).shuffle(cells)
    rows = []
    for j, (pdb, seed, cap) in enumerate(cells):
        directory = args.work / pdb
        if json.loads((directory / "screening.json").read_text())["status"] != "ready":
            raise ValueError("Ineligible " + pdb)
        meta = json.loads((directory / "preparation.json").read_text())
        order = list(arms)
        if args.mode == "persistent" and cap in (2000, 9000):
            order = ["lock"]
        cut = j % len(order)
        for arm in order[cut:] + order[:cut]:
            key = f"{pdb}_{seed}_{cap}_{arm}"
            dest = raw / (key + ".pdbqt")
            box = [
                x
                for d, center in zip("xyz", meta["center"])
                for x in ("--center_" + d, str(center), "--size_" + d, "26")
            ]
            opts = [
                "--receptor",
                str(directory / "receptor.pdbqt"),
                "--ligand",
                str(directory / "ligand.pdbqt"),
                *box,
                "--cpu",
                "1",
                "--seed",
                str(seed),
                "--max_evals",
                str(cap),
                "--num_modes",
                "5",
                "--verbosity",
                "0",
                "--exhaustiveness",
                "1",
                "--energy_range",
                "100",
                "--out",
                str(dest),
            ]
            start = time.perf_counter_ns()
            log = subprocess.check_output(
                [
                    arms[arm],
                    *opts,
                    *([] if arm == "upstream" else ["--search_counters"]),
                ],
                stderr=subprocess.STDOUT,
                text=True,
            )
            elapsed = time.perf_counter_ns() - start
            (raw / (key + ".log")).write_text(log)
            if arm == "upstream":
                replay = raw / (key + ".verify")
                log = subprocess.check_output(
                    [args.upstream_counts, *opts[:-1], replay],
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                (raw / (key + ".counts")).write_text(log)
                if sha256(replay) != sha256(dest):
                    raise ValueError("Upstream changed under counting replay: " + key)
                replay.unlink()
            counters = {
                k: int(v)
                for k, v in re.findall(
                    r"(hunt_calls|hunt_evals|second_calls|second_evals|explicit_calls|explicit_evals)=(\d+)",
                    log,
                )
            }
            if len(counters) != 6:
                raise ValueError(key + ": missing counters")
            rows.append(
                dict(
                    key=key,
                    pdb=pdb,
                    arm=arm,
                    seed=seed,
                    cap=cap,
                    elapsed_ns=elapsed,
                    **counters,
                    joint=sum(v for k, v in counters.items() if k.endswith("evals")),
                    modes=dest.read_text().count("ENDMDL"),
                    sha256=sha256(dest),
                )
            )
            write(output / "runs.csv", rows)
    (output / "binary_hashes.json").write_text(
        json.dumps(
            {
                k: sha256(v)
                for k, v in {**arms, "upstream_counts": args.upstream_counts}.items()
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for field in (
        "work",
        "output",
        "upstream",
        "upstream_counts",
        "parent",
        "reset",
        "ccd",
    ):
        p.add_argument("--" + field.replace("_", "-"), type=Path, required=True)
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    p.add_argument("--lock", type=Path)
    p.add_argument(
        "--mode", choices=["paired", "lower", "repeat", "persistent"], default="paired"
    )
    p.add_argument("--ids", nargs="+")
    run(p.parse_args())
