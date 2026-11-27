"""Diagnostic replay of final MC pool, promotion filter and explicit refinement."""

import argparse
import csv
import json
import os
import sys
from pathlib import Path

from backtracking_panel import BENCH, CASES, paired

sys.path.insert(0, str(BENCH / "pair_gradient/holdout"))
from run_holdout import stage_data


def run(work, binary, cpu):
    os.sched_setaffinity(0, {cpu})
    w = work / "trace"
    w.mkdir(exist_ok=True)
    (w / "runs").mkdir(exist_ok=True)
    meta = {q: json.loads((work / q / "preparation.json").read_text()) for q in CASES}
    for q in CASES:
        if not (w / q).exists():
            (w / q).symlink_to(work / q, target_is_directory=True)
    with (work / "cells.csv").open() as f:
        rows = [r for r in csv.DictReader(f) if r["arm"] != "upstream"]
    data = []
    counts = []
    for r in rows:
        arm = r["arm"]
        q = r["query"]
        rx = r["receiver"]
        os.environ["FASTDOCK_RESEARCH_BACKTRACK"] = "native" if arm == "parent" else arm
        rr, _ = paired.execute(
            w, meta, q, rx, int(r["seed"]), int(r["cap"]), arm, binary, False
        )
        key = rr["key"]
        out = key + ".pdbqt"
        if (w / "runs" / out).read_bytes() != (work / "runs" / out).read_bytes():
            raise ValueError("Passive stage trace changed poses in " + key)
        raw = stage_data(
            (w / "runs" / (key + ".log")).read_text(), key, meta[q], meta[rx]
        )
        data.extend(raw)
        stages = {
            s: [p for p in raw if p["stage"] == s]
            for s in ("search", "before", "after")
        }
        if len(stages["before"]) != len(stages["after"]) or len(
            stages["before"]
        ) != int(r["explicit_inputs"]):
            raise ValueError("Refinement input count mismatch in " + key)
        counts.append(
            {
                "key": key,
                "arm": arm,
                "target": r["target"],
                **{s: min(p["rmsd_A"] for p in poses) for s, poses in stages.items()},
                "returned": r["oracle_A"],
                "top1": r["top1_A"],
            }
        )
    paired.write(work / "stages.csv", data)
    paired.write(work / "stage_oracles.csv", counts)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--cpu", required=True, type=int)
    a = p.parse_args()
    run(a.work.resolve(), a.binary.resolve(), a.cpu)
