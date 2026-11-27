"""Replay MC saving, input selection and explicit refinement without changing any pose."""

import argparse
import csv
import json
import os
import sys
from pathlib import Path

from panel import STUDIES, configure_arm, design, digest, graph_and_receptor_tools


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("observer", type=Path)
    p.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3])
    p.add_argument("--cpu", type=int, required=True)
    p.add_argument("--study", choices=STUDIES, default="tangent")
    a = p.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    cases, _, _ = design(a.study)
    panel, _ = graph_and_receptor_tools(a.repo.resolve(), cases)
    sys.path.insert(0, str(a.repo / "benchmarks/pair_gradient/holdout"))
    from run_holdout import stage_data

    work = a.work.resolve()
    meta = {q: json.loads((work / q / "preparation.json").read_text()) for q in cases}
    trace = work / "trace"
    (trace / "runs").mkdir(parents=True, exist_ok=True)
    for q in cases:
        if not (trace / q).exists():
            (trace / q).symlink_to(work / q, target_is_directory=True)
    with (work / "cells.csv").open() as handle:
        cells = list(csv.DictReader(handle))
    stages, totals = [], []
    for row in cells:
        arm = row["arm"]
        if arm == "upstream":
            continue
        configure_arm(arm)
        q, r, key = row["query"], row["receiver"], row["key"]
        replay, _ = panel.paired.execute(
            trace, meta, q, r, int(row["seed"]), int(row["cap"]), arm, a.observer, False
        )
        if digest(trace / "runs" / (key + ".pdbqt")) != row["pdbqt_sha256"]:
            raise ValueError("Observer changed output: " + key)
        log = (trace / "runs" / (key + ".log")).read_text()
        if any(str(replay[k]) != row[k] for k in panel.paired.counts(log)):
            raise ValueError("Observer changed work counts: " + key)
        entries = stage_data(log, key, meta[q], meta[r])
        groups = {
            name: [v for v in entries if v["stage"] == name]
            for name in ("search", "before", "after")
        }
        if not (
            len(groups["after"]) == len(groups["before"]) == int(row["explicit_inputs"])
        ):
            raise ValueError("Mismatched explicit refinement input count: " + key)
        stages.extend(entries)
        totals.append(
            dict(
                key=key,
                arm=arm,
                target=row["target"],
                **{s: min(p["rmsd_A"] for p in group) for s, group in groups.items()},
                returned=row["oracle_A"],
                top1=row["top1_A"],
            )
        )
    panel.paired.write(work / "stages.csv", stages)
    panel.paired.write(work / "stage_oracles.csv", totals)


if __name__ == "__main__":
    main()
