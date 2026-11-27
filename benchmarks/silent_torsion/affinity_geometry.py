"""Measure receptor overlaps and score–size dependence at the actual Gal3 first poses."""

import argparse
import csv
from pathlib import Path
import re
import sys

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH / "rigid_flip"))
from holdout import paired


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def geometry(work, inputs):
    description = {
        r["compound"]: r
        for r in read(BENCH / "pair_gradient/galectin_external/manifest.csv")
    }
    trees = {}
    for receptor in ("1KJL", "3T1L"):
        atoms = paired.atoms((inputs / (receptor + ".pdbqt")).read_text().splitlines())
        trees[receptor] = (
            cKDTree(np.array([xyz for xyz, r in atoms.values()])),
            np.array([r for xyz, r in atoms.values()]),
        )
    rows = []
    for cell in read(work / "cells.csv"):
        receptor, lig, arm = [cell[k] for k in ("receptor", "compound", "arm")]
        file = work / "runs" / f"{receptor}_{lig}_{arm}.pdbqt"
        pose = file.read_text().split("ENDMDL", 1)[0]
        atoms = paired.atoms(pose.splitlines())
        xyz = np.array([x for x, r in atoms.values()])
        radii = np.array([r for x, r in atoms.values()])
        tree, rr = trees[receptor]
        overlaps = 0
        for k, close in enumerate(tree.query_ball_point(xyz, 3.2)):
            if close:
                d = np.linalg.norm(xyz[k] - tree.data[close], axis=1)
                overlaps += int(np.sum(d < 0.75 * (radii[k] + rr[close])))
        terms = dict(
            re.findall(r"^REMARK (INTER|INTRA|UNBOUND):\s*([-.\d]+)", pose, re.M)
        )
        rows.append(
            dict(
                arm=arm,
                receiver=receptor,
                compound=lig,
                heavy=int(description[lig]["heavy"]),
                score=float(cell["score"]),
                inter=float(terms["INTER"]),
                intra=float(terms["INTRA"]),
                unbound=float(terms["UNBOUND"]),
                clashes=overlaps,
                min_distance_A=float(tree.query(xyz)[0].min()),
            )
        )
    paired.write(work / "top_geometry.csv", rows)
    summary = []
    for arm in sorted({r["arm"] for r in rows}):
        group = [r for r in rows if r["arm"] == arm]
        averages = [
            np.mean([r["score"] for r in group if r["compound"] == lig])
            for lig in description
        ]
        sizes = [int(r["heavy"]) for r in description.values()]
        kd = [float(r["kd_uM"]) for r in description.values()]
        summary.append(
            dict(
                arm=arm,
                poses=len(group),
                clashing_poses=sum(r["clashes"] > 0 for r in group),
                clash_pairs=sum(r["clashes"] for r in group),
                score_vs_heavy_spearman=spearmanr(averages, sizes).statistic,
                kd_vs_heavy_spearman=spearmanr(kd, sizes).statistic,
                max_intra_unbound_delta=max(
                    abs(r["intra"] - r["unbound"]) for r in group
                ),
            )
        )
    paired.write(work / "geometry_summary.csv", summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("inputs", type=Path)
    a = p.parse_args()
    geometry(a.work, a.inputs)
