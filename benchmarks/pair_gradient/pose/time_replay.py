"""Repeat only timings of the frozen held-out jobs; require byte-identical coordinates."""

import argparse
import csv
import json
import random
import resource
import subprocess
import tempfile
import time
from pathlib import Path

from paired import write


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    for a in ("parent", "linear", "hermite", "shift"):
        ap.add_argument("--" + a, required=True, type=Path)
    args = ap.parse_args()
    w = args.work.resolve()
    with (w / "cells.csv").open() as f:
        cells = [
            c
            for c in csv.DictReader(f)
            if c["arm"] == "parent" and c["seed"] == "98251" and c["cap"] == "12000"
        ]
    with (w / "matched_cells.csv").open() as f:
        cells.extend(csv.DictReader(f))
    results = []
    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "pose.pdbqt"
        for repeat in range(1, 4):
            order = sorted(cells, key=lambda c: c["key"])
            random.Random(20261127 + repeat).shuffle(order)
            for c in order:
                query, receiver, arm = c["query"], c["receiver"], c["arm"]
                meta = json.loads((w / receiver / "preparation.json").read_text())
                cmd = [
                    str(getattr(args, arm).resolve()),
                    "--receptor",
                    str(w / receiver / "receptor.pdbqt"),
                    "--ligand",
                    str(w / query / "ligand.pdbqt"),
                    "--cpu",
                    "1",
                    "--exhaustiveness",
                    "1",
                    "--max_evals",
                    c["cap"],
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--seed",
                    c["seed"],
                    "--out",
                    str(out),
                    "--search_counters",
                ]
                for axis, val in zip("xyz", meta["center"]):
                    cmd += ["--center_" + axis, str(val), "--size_" + axis, "24"]
                before_cpu = resource.getrusage(resource.RUSAGE_CHILDREN)
                before = time.perf_counter()
                subprocess.run(
                    cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE
                )
                seconds = time.perf_counter() - before
                after_cpu = resource.getrusage(resource.RUSAGE_CHILDREN)
                cpu = (
                    after_cpu.ru_utime
                    - before_cpu.ru_utime
                    + after_cpu.ru_stime
                    - before_cpu.ru_stime
                )
                if (
                    out.read_bytes()
                    != (w / "runs" / (c["key"] + ".pdbqt")).read_bytes()
                ):
                    raise ValueError(
                        "timing-only replay failed coordinate parity " + c["key"]
                    )
                results.append(
                    {
                        "key": c["key"],
                        "repeat": repeat,
                        "seconds": seconds,
                        "cpu_seconds": cpu,
                        "byte_identical": True,
                    }
                )
            print("Timing repeat", repeat, flush=True)
    write(w / "time_replay.csv", results)


if __name__ == "__main__":
    main()
