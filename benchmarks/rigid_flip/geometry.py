"""Independently check heavy-atom RMSD and quantify centroid, orientation and shape errors."""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from holdout import BENCH, CASES, paired

sys.path.insert(0, str(BENCH / "bond_distance/pose"))
from audit import independent_rmsds

sys.path.insert(0, str(BENCH / "paired_crossdocking"))
from crystal_start import run as crystal_starts


def errors(query, rec, rows, path):
    rot, trans, _, _ = paired.base.alignment(query, rec)
    refs = np.asarray(query["mapped"]) @ rot + trans
    checked = independent_rmsds(query, rec, path)
    blocks = path.read_text().split("ENDMDL")[:-1]
    if not (len(blocks) == len(checked) == len(rows)):
        raise ValueError(path)
    out = []
    for i, (block, reference_rmsd, row) in enumerate(zip(blocks, checked, rows), 1):
        atoms = paired.atoms(block.splitlines())
        xyz = np.array([atoms[j][0] for j in query["serial"]])
        distances = np.square(refs - xyz).sum(axis=2).mean(axis=1)
        ref = refs[np.argmin(distances)]
        center_sq = np.sum((xyz.mean(0) - ref.mean(0)) ** 2)
        x, y = xyz - xyz.mean(0), ref - ref.mean(0)
        u, _, vt = np.linalg.svd(x.T @ y)
        u[:, -1] *= np.sign(np.linalg.det(u @ vt))
        shape_sq = np.square(x @ u @ vt - y).sum(axis=1).mean()
        orient_sq = np.square(x - y).sum(axis=1).mean() - shape_sq
        if (
            abs(float(row["rmsd_A"]) - reference_rmsd) > 1e-7
            or abs(center_sq + shape_sq + orient_sq - reference_rmsd**2) > 1e-7
        ):
            raise ValueError((path, i))
        out.append(
            {
                "key": row["key"],
                "mode": i,
                "independent_rmsd_A": reference_rmsd,
                "centroid_A": np.sqrt(center_sq),
                "internal_A": np.sqrt(shape_sq),
                "orientation_A": np.sqrt(max(0, orient_sq)),
            }
        )
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument(
        "--crystal", type=Path, help="Pro binary for offline starting-pose diagnostic"
    )
    a = ap.parse_args()
    w = a.work.resolve()
    meta = {q: json.loads((w / q / "preparation.json").read_text()) for q in CASES}
    with (w / "modes.csv").open() as f:
        by_key = {}
        for row in csv.DictReader(f):
            by_key.setdefault(row["key"], []).append(row)
    pairs = set()
    result = []
    with (w / "cells.csv").open() as f:
        for row in csv.DictReader(f):
            key = row["key"]
            q, r = row["query"], row["receiver"]
            pairs.add((q, r))
            result.extend(
                errors(meta[q], meta[r], by_key[key], w / "runs" / (key + ".pdbqt"))
            )
    paired.write(w / "geometry.csv", result)
    if a.crystal:
        import sources

        sources.CASES = CASES
        sources.reference_geometry(w)
        (w / "metadata.json").write_text(json.dumps(meta))
        crystal_starts(w, a.crystal.resolve(), pairs, compare_mappings=True)
