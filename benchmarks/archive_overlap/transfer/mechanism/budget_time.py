"""Counterbalance three single-CPU, uninstrumented baseline/reduced-cap timing repeats."""

import argparse
import csv
import json
import os
from pathlib import Path
import random

from transfer_stages import arguments, run_panel

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("panel", "reduced", "output", "parent", "candidate"):
        p.add_argument(name, type=Path)
    p.add_argument("--cpu", type=int, required=True)
    a = p.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    source = a.panel.resolve()
    dest = a.output.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    rows = [
        x
        for x in csv.DictReader((a.reduced / "cells.csv").open())
        if x["arm"] == "candidate"
    ]
    meta = {
        r["query"]: json.loads((source / r["query"] / "preparation.json").read_text())
        for r in rows
    }
    results = []
    for repeat in range(3):
        order = rows.copy()
        random.Random(2287 + repeat).shuffle(order)
        for i, row in enumerate(order):
            q, r, seed = row["query"], row["receiver"], row["seed"]
            cases = [
                ("parent", 12000, source),
                ("candidate", 11400, a.reduced.resolve()),
            ]
            if (i + repeat) % 2:
                cases.reverse()
            for arm, cap, reference in cases:
                key = f"{q}_{r}_{seed}_{arm}"
                common = arguments(source, meta, q, r, seed, cap)
                output, log = dest / (key + ".pdbqt"), dest / (key + ".log")
                wall, cpu, rss = run_panel.cli(
                    getattr(a, arm).resolve(), common, output, log
                )
                if (
                    output.read_bytes()
                    != (reference / "runs" / output.name).read_bytes()
                ):
                    raise ValueError(key + " differs")
                results.append(
                    dict(
                        repeat=repeat + 1,
                        key=key,
                        arm=arm,
                        cap=cap,
                        wall_s=wall,
                        cpu_s=cpu,
                        maxrss_KiB=rss,
                    )
                )
        run_panel.paired.write(dest / "paired_cost.csv", results)
