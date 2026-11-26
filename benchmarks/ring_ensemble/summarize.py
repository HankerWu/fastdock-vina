"""Tabulate per-target structure outcomes and input-to-crystal local ring errors."""

import argparse
import csv
import json
import sys
from pathlib import Path
from statistics import median

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from ring_conformers import ring_distance
from run import write_csv


def summarize(work, tables):
    with (tables / "cases.csv").open() as f:
        rows = list(csv.DictReader(f))
    counts = []
    for group in [
        "all",
        *sorted({r["target"] for r in rows}),
        *sorted({r["kind"] for r in rows}),
    ]:
        for arm in ["first", "repeat", "ring", "sham"]:
            cell = [
                r
                for r in rows
                if r["arm"] == arm and group in ("all", r["target"], r["kind"])
            ]
            counts.append(
                {
                    "group": group,
                    "arm": arm,
                    "cells": len(cell),
                    "top_hits": sum(int(r["top_hit"]) for r in cell),
                    "pool_hits": sum(int(r["pool_hit"]) for r in cell),
                    "naive_top_hits": sum(float(r["naive_rmsd"]) <= 2 for r in cell),
                    "gas_offset_top_hits": sum(
                        float(r["gas_offset_rmsd"]) <= 2 for r in cell
                    ),
                    "median_top_rmsd": median(float(r["top_rmsd"]) for r in cell),
                }
            )
    write_csv(tables / "summary.csv", counts)
    rows = []
    for p in sorted((work / "prepared").glob("*/panel.json")):
        data = json.loads(p.read_text())
        states = data["states"]
        ref = np.array(data["receptors"][data["crystal"]]["reference"])
        row = {
            "target": data["target"],
            "ccd": data["ccd"],
            "ring_distinct": states["alternate_distinct"],
            "sham_distinct": states["sham_distinct"],
            "separation_A": states["ring_distance_A"],
            "alternate_gas_delta": np.diff(states["energies_MMFF94s"])[0],
            "first_gas_delta": states["first_energy_MMFF94s"]
            - states["energies_MMFF94s"][0],
        }
        for state in ["first", "base", "alternate", "sham"]:
            ids = {}
            pos = {}
            for line in (p.parent / (state + ".pdbqt")).read_text().splitlines():
                if line.startswith("REMARK INDEX MAP"):
                    v = list(map(int, line.split()[3:]))
                    ids.update(zip(v[::2], v[1::2], strict=True))
                elif line.startswith(("ATOM", "HETATM")):
                    pos[int(line[6:11])] = [
                        float(line[i : i + 8]) for i in (30, 38, 46)
                    ]
            xyz = np.array([pos[ids[i + 1]] for i in range(len(ref))])
            row[state] = ring_distance(
                ref, xyz, states["fragments"], data["permutations"]
            )
        rows.append(row)
    write_csv(tables / "ring_geometry.csv", rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("work", type=Path)
    parser.add_argument("tables", type=Path)
    args = parser.parse_args()
    summarize(args.work, args.tables)
