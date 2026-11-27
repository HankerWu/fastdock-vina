"""Paired replay of the lowest-energy and longest-visit selectors on common trajectories."""

import argparse
import csv
import hashlib
import os
import random
import re
import subprocess
import sys
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "coverage"))
from probe import PANEL, command, counters, mode_distances, references, table


def jobs(args):
    historical = list(csv.DictReader((PANEL / "cells.csv").open()))
    for r in historical:
        yield args.archived, r, "development"
    seen = set()
    for r in historical:
        pair = r["query"], r["receiver"]
        if pair in seen:
            continue
        seen.add(pair)
        yield args.archived, dict(r, seed="411119"), "new_seed"
    for family, ids in (("NUDT1", ("4N1U", "4C9X")), ("CHIT1", ("5NRA", "6ZE8"))):
        for q in ids:
            for r in ids:
                yield (
                    args.structures,
                    {
                        "target": family,
                        "query": q,
                        "receiver": r,
                        "seed": "411119",
                        "kind": "cognate" if q == r else "cross",
                    },
                    "new_seed",
                )


def run(args):
    os.sched_setaffinity(0, {args.cpu})
    args.out.mkdir(exist_ok=True, parents=True)
    raw = args.out / "runs"
    raw.mkdir(exist_ok=True)
    gen = random.Random(411113)
    experiments = list(jobs(args))
    gen.shuffle(experiments)
    rows = []
    binaries = {
        "parent": args.parent,
        "energy": args.energy,
        "recurrence": args.recurrence,
    }
    for inputs, case, split in experiments:
        serial, refs = references(inputs, case)
        arms = list(binaries)
        gen.shuffle(arms)
        observed = []
        pending = {}
        for arm in arms:
            key = (
                "_".join(case[k] for k in ("target", "query", "receiver", "seed"))
                + "_"
                + arm
            )
            out = raw / (key + ".pdbqt")
            start = time.perf_counter_ns()
            p = subprocess.run(
                command(inputs, case, binaries[arm], out),
                capture_output=True,
                check=True,
            )
            ns = time.perf_counter_ns() - start
            (raw / (key + ".log")).write_bytes(p.stdout + p.stderr)
            ct = counters(p.stdout)
            assert len(ct) == 4
            costs = dict(re.findall(rb"(\w+)=(\d+)", b" ".join(ct)))
            deferred = re.search(rb"Deferred hunt: ([^\r\n]+)", p.stdout)
            tags = dict(re.findall(rb"(\w+)=([^ ]+)", deferred[1])) if deferred else {}
            assert bool(tags) == (arm != "parent")
            scores = mode_distances(out.read_bytes(), serial, refs)
            row = {k: case[k] for k in ("target", "kind", "query", "receiver", "seed")}
            row.update(
                split=split,
                arm=arm,
                key=key,
                ns=ns,
                joint=sum(
                    int(costs[k])
                    for k in (b"second_evals", b"hunt_evals", b"explicit_evals")
                ),
                local_calls=sum(
                    int(costs[k])
                    for k in (b"second_calls", b"hunt_calls", b"explicit_calls")
                ),
                slots=int(costs[b"explicit_inputs"]),
                replaced=int(tags.get(b"replaced", 0)),
                visit=int(tags.get(b"run_length", 0)),
                top1=scores[0],
                top5=min(scores[:5]),
                oracle=min(scores),
                sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
            )
            rows.append(row)
            observed.append((ct[:3], row["slots"]))
            pending[arm] = row
        if any(x != observed[0] for x in observed):
            raise ValueError("MC or explicit input quota diverged")
        for arm in ("energy", "recurrence"):
            if (
                not pending[arm]["replaced"]
                and pending[arm]["sha256"] != pending["parent"]["sha256"]
            ):
                raise ValueError("Unexpected changed output")
    table(args.out / "cells.csv", rows)
    with tarfile.open(args.out / "runs.tar.xz", "w:xz") as t:
        for f in sorted(raw.iterdir()):
            t.add(f, arcname=f.name)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("archived", "structures", "out", "parent", "energy", "recurrence"):
        p.add_argument(k, type=Path)
    p.add_argument("--cpu", type=int, required=True)
    run(p.parse_args())
