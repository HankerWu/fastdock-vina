"""Replay equal-input searches in randomized order and require byte-identical output."""

import argparse
import csv
import json
import os
import random
import resource
import subprocess
import tempfile
import time
from pathlib import Path
from inputs import reference


def replay(work, binaries):
    with (work / "cells.csv").open() as f:
        cells = [c for c in csv.DictReader(f) if c["arm"] in binaries]
    rows = []
    with tempfile.TemporaryDirectory() as temp:
        out = Path(temp) / "out.pdbqt"
        for n in range(1, 4):
            order = sorted(cells, key=lambda c: c["key"])
            random.Random(214087 + n).shuffle(order)
            for c in order:
                q, r, arm = c["query"], c["receiver"], c["arm"]
                cmd = [
                    binaries[arm],
                    "--receptor",
                    work / r / "receptor.pdbqt",
                    "--ligand",
                    work / q / "ligand.pdbqt",
                    "--out",
                    out,
                    "--exhaustiveness",
                    "1",
                    "--cpu",
                    "1",
                    "--energy_range",
                    "100",
                    "--seed",
                    c["seed"],
                    "--num_modes",
                    "9",
                    "--max_evals",
                    c["cap"],
                    "--search_counters",
                ]
                for axis, x in zip(
                    "xyz",
                    json.loads((work / r / "preparation.json").read_text())["center"],
                ):
                    cmd.extend(["--center_" + axis, str(x), "--size_" + axis, "24"])
                before = resource.getrusage(resource.RUSAGE_CHILDREN)
                start = time.perf_counter()
                subprocess.run(
                    cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=True
                )
                seconds = time.perf_counter() - start
                after = resource.getrusage(resource.RUSAGE_CHILDREN)
                cpu = (after.ru_utime - before.ru_utime) + (
                    after.ru_stime - before.ru_stime
                )
                if (
                    out.read_bytes()
                    != (work / "runs" / (c["key"] + ".pdbqt")).read_bytes()
                ):
                    raise ValueError(c["key"])
                rows.append(
                    dict(
                        key=c["key"],
                        repeat=n,
                        seconds=seconds,
                        cpu_seconds=cpu,
                        output_identical=True,
                    )
                )
            print(n, flush=True)
    reference.write(work / "timing_replays.csv", rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument("--cpu", type=int, default=8)
    arms = ("parent", "gated", "unconditional")
    for k in arms:
        ap.add_argument("--" + k, type=Path, required=True)
    a = ap.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    replay(a.work, {k: getattr(a, k).resolve() for k in arms})
