"""Test deployable selection on exactly the same six-family trajectories as the coverage audit."""

import argparse
import csv
import hashlib
import os
import random
import re
import subprocess
import tarfile
import time
from pathlib import Path

from probe import command, counters, mode_distances, references, table


def main(a):
    os.sched_setaffinity(0, {a.cpu})
    cases = list(csv.DictReader((a.results / "cells.csv").open()))
    gen = random.Random(66829)
    gen.shuffle(cases)
    files = a.work / "blind"
    files.mkdir(exist_ok=True)
    rows = []
    for case in cases:
        base = a.work / "runs" / case["key"]
        nativecounts = counters(base.with_suffix(".log").read_bytes())
        serial, refs = references(a.work, case)
        arms = [("parent", a.parent), ("gated", a.gated), ("ungated", a.ungated)]
        gen.shuffle(arms)
        for arm, binary in arms:
            dest = files / (case["key"] + "_" + arm + ".pdbqt")
            tic = time.perf_counter_ns()
            call = subprocess.run(
                command(a.work, case, binary, dest), check=True, capture_output=True
            )
            elapsed = time.perf_counter_ns() - tic
            data = dest.read_bytes()
            current = counters(call.stdout)
            if len(current) != 4 or current[:3] != nativecounts[:3]:
                raise ValueError("Change to global search or pool merge")
            fields = dict(re.findall(rb"(\w+)=(\d+)", b" ".join(current)))
            if int(fields[b"explicit_inputs"]) != int(
                re.search(rb"explicit_inputs=(\d+)", nativecounts[-1])[1]
            ):
                raise ValueError("Different number of refinement slots")
            observe = re.search(rb"Deferred hunt: (.*)$", call.stdout, re.MULTILINE)
            tags = (
                {}
                if observe is None
                else dict(re.findall(rb"(\w+)=([^ ]+)", observe[1]))
            )
            if bool(tags) != (arm != "parent"):
                raise ValueError("Missing audit flag")
            sha = hashlib.sha256(data).hexdigest()
            replaced = int(tags.get(b"replaced", 0))
            if not replaced and sha != case["native_sha256"]:
                raise ValueError("Unchanged slot allocation altered output")
            scores = mode_distances(data, serial, refs)
            row = {
                k: case[k]
                for k in ("target", "kind", "query", "receiver", "seed", "key")
            }
            row.update(
                arm=arm,
                ns=elapsed,
                joint=sum(
                    int(fields[k])
                    for k in (b"second_evals", b"hunt_evals", b"explicit_evals")
                ),
                local_calls=sum(
                    int(fields[k])
                    for k in (b"second_calls", b"hunt_calls", b"explicit_calls")
                ),
                explicit=int(fields[b"explicit_evals"]),
                slots=int(fields[b"explicit_inputs"]),
                eligible=int(tags.get(b"eligible", 0)),
                replaced=replaced,
                slot=int(tags.get(b"slot", 0)),
                top1=scores[0],
                top5=min(scores[:5]),
                oracle=min(scores),
                sha256=sha,
            )
            rows.append(row)
            dest.with_suffix(".log").write_bytes(call.stdout + call.stderr)
            if sha == case["native_sha256"]:
                dest.unlink()
    table(a.results / "blind.csv", rows)
    with tarfile.open(a.results / "blind.tar.xz", "w:xz") as t:
        for f in sorted(files.iterdir()):
            t.add(f, arcname=f.name)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "results", "parent", "gated", "ungated"):
        p.add_argument(name, type=Path)
    p.add_argument("--cpu", type=int, required=True)
    main(p.parse_args())
