"""Verify both real-data mappings, compute paired changes, and compare initial affected minimizations."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from replay import evaluate, write_csv


def trace_events(path):
    lines = path.read_text().splitlines()
    if len(lines) % 3:
        raise ValueError(f"Incomplete event in {path}")
    events = []
    for i in range(0, len(lines), 3):
        start, next_step, end = (line.split(",") for line in lines[i : i + 3])
        if (start[0], next_step[0], end[0]) != ("first", "next", "finish") or len(
            {start[1], next_step[1], end[1]}
        ) != 1:
            raise ValueError(f"Inconsistent event in {path}")
        events.append((start, next_step, end))
    return events


def check(work):
    meta = json.loads((work / "metadata.json").read_text())
    with (work / "cells.csv").open() as stream:
        cells = list(csv.DictReader(stream))
    with (work / "poses.csv").open() as stream:
        modes = {
            (r["case"], r["cap"], r["seed"], r["arm"], int(r["rank"])): float(
                r["rmsd_A"]
            )
            for r in csv.DictReader(stream)
        }
    rows, seen, max_error = [], {}, 0.0
    for r in cells:
        key = (r["case"], r["cap"], r["seed"], r["arm"])
        tag = "_".join(key)
        file = work / "runs" / (tag + ".pdbqt")
        if hashlib.sha256(file.read_bytes()).hexdigest() != r["pool_sha256"]:
            raise ValueError(f"Coordinate hash mismatch in {file}")
        # Second mapping path: use archived CCD coordinates/serials instead of
        # inferring the deposition's stereochemical graph from the PDB block again.
        case = meta[r["case"]]
        mapped = [
            evaluate(file, case["serial"], np.asarray(ref)) for ref in case["mapped"]
        ]
        for equivalent in zip(*mapped):
            rank = equivalent[0]["rank"]
            rmsd = min(m["rmsd_A"] for m in equivalent)
            max_error = max(max_error, abs(modes.pop((*key, rank)) - rmsd))
        seen[key] = r
    if modes or max_error > 1e-6:
        raise ValueError(
            f"Unmatched modes or discrepant mappings: {len(modes)}, {max_error}"
        )

    events_count = {"parent": 0, "fixed": 0}
    for (case, cap, seed, arm), fixed in seen.items():
        if arm != "fixed":
            continue
        key = (case, cap, seed)
        parent = seen[(*key, "parent")]
        row = dict(case=case, cap=cap, seed=seed)
        for field in ("top1_A", "best5_A", "oracle_A", "regret_A"):
            row["delta_" + field] = float(fixed[field]) - float(parent[field])
        events = {}
        for mode in events_count:
            trace = work / "traces" / ("_".join((*key, mode)) + ".trace")
            events[mode] = trace_events(trace)
            events_count[mode] += len(events[mode])
            row[mode + "_events"] = len(events[mode])
        extra = {
            "first_sy": "",
            "parent_pg": "",
            "fixed_pg": "",
            "first_accepted_energy": "",
            "parent_returned": "",
            "fixed_returned": "",
            "parent_FG": "",
            "fixed_FG": "",
        }
        if events["parent"] or events["fixed"]:
            old = events["parent"][0]
            new = events["fixed"][0]
            if old[0] != new[0]:
                raise ValueError(f"First affected start differs in {key}")
            extra.update(
                first_sy=float(old[0][6]),
                parent_pg=float(old[1][2]),
                fixed_pg=float(new[1][2]),
                first_accepted_energy=float(old[0][4]),
                parent_returned=float(old[2][3]),
                fixed_returned=float(new[2][3]),
                parent_FG=int(old[2][2]),
                fixed_FG=int(new[2][2]),
            )
        rows.append({**row, **extra})
    write_csv(work / "paired.csv", rows)
    stats = dict(
        cells=len(cells),
        matched_cells=len(rows),
        poses=sum(int(r["modes"]) for r in cells),
        max_mapping_discrepancy_A=max_error,
        total_affected_minimizations=events_count,
    )
    (work / "verification.json").write_text(json.dumps(stats, indent=2) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    check(ap.parse_args().work)
