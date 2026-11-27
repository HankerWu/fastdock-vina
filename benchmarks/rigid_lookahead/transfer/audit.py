"""Recompute every output RMSD from coordinates, and separate rigid and internal shape errors."""

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from prepare_transfer import DIRECTIONS, SYSTEMS


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def read_atoms(block, serial):
    points = {int(line[6:11]): [float(line[j:j + 8]) for j in (30, 38, 46)]
              for line in block.splitlines() if line.startswith(("ATOM", "HETATM"))}
    return np.asarray([points[i] for i in serial])


def reference(query, receiver):
    common = sorted(set(query["ca"]).intersection(receiver["ca"]))
    x, y = (np.asarray([m["ca"][k] for k in common]) for m in (query, receiver))
    keep = (np.linalg.norm(x - query["center"], axis=1) < 14) & (
        np.linalg.norm(y - receiver["center"], axis=1) < 14)
    x, y = x[keep], y[keep]
    fit, rss = Rotation.align_vectors(y - y.mean(0), x - x.mean(0))
    refs = np.asarray([fit.apply(np.asarray(a) - x.mean(0)) + y.mean(0)
                       for a in query["mapped"]])
    return refs, float(rss / np.sqrt(len(x))), len(x)


def check(base, baseline, rigid):
    meta = {c: json.loads((base / c / "preparation.json").read_text()) for c in SYSTEMS}
    pairs = {(r, q) for r, q in DIRECTIONS if meta[r]["receptor_status"] == "ready"}
    refs = {key: reference(meta[key[1]], meta[key[0]]) for key in pairs}
    shapes = {}
    for code in sorted({q for _, q in pairs}):
        xyz = read_atoms((base / code / "ligand.pdbqt").read_text(), meta[code]["serial"])
        fixed = xyz - xyz.mean(0)
        shapes[code] = min(float(Rotation.align_vectors(
            np.asarray(ref) - np.mean(ref, axis=0), fixed)[1] / np.sqrt(len(xyz)))
            for ref in meta[code]["mapped"])

    errors = []
    def audit_model(receiver, query, block, recorded):
        xyz = read_atoms(block, meta[query]["serial"])
        measured = np.min(np.sqrt(np.sum((refs[receiver, query][0] - xyz) ** 2, axis=(1, 2)) / len(xyz)))
        errors.append(abs(float(recorded) - measured))

    base_rows = table(baseline / "baseline.csv")
    base_modes = table(baseline / "baseline_poses.csv")
    if len(base_rows) != len(pairs) * 8:
        raise ValueError("Missing baseline cells")
    for row in base_rows:
        key = row["key"]
        data = (baseline / "raw" / (key + ".pdbqt")).read_bytes()
        if hashlib.sha256(data).hexdigest() != row["output_sha256"]:
            raise ValueError("Changed baseline output: " + key)
        models = [s for s in data.decode().split("ENDMDL") if "MODEL" in s]
        items = [a for a in base_modes if a["key"] == key]
        if len(items) != len(models) or len(items) != int(row["output_poses"]):
            raise ValueError("Missing baseline pose: " + key)
        for m, mode in zip(models, items):
            audit_model(row["receiver"], row["query"], m, mode["rmsd"])

    poses = table(rigid / "poses.csv")
    cells = table(rigid / "cells.csv")
    ranks = table(rigid / "ranked.csv")
    keys = ("receiver", "query", "seed", "arm")
    def key(row):
        return tuple(row[k] for k in keys)

    for a in poses:
        r, q, seed, arm = key(a)
        start = rigid / "runs" / "_".join(key(a)) / (a["slot"] + ".pdbqt")
        audit_model(r, q, start.read_text(), a["before_rmsd"])
        audit_model(r, q, Path(str(start) + ".refined").read_text(), a["after_rmsd"])

    for cell in cells:
        ranked = sorted((r for r in ranks if key(cell) == key(r)), key=lambda r: int(r["rank"]))
        vals = [float(r["after_rmsd"]) for r in ranked]
        if (len(vals) != int(cell["reported"]) or
            [vals[0], min(vals[:5]), min(vals)] != [float(cell[k]) for k in ("top", "top5", "oracle")]):
            raise ValueError("Incorrect output summary: " + str(key(cell)))
    conditions = set(itertools.product(pairs, ("43103", "43117"), ("raw", "sham", "gradient")))
    if {(key(c)[:2], c["seed"], c["arm"]) for c in cells} != conditions:
        raise ValueError("Incorrect transfer cells")
    if max(errors) > 1e-8:
        raise ValueError("RMSD replay discrepancy: " + str(max(errors)))
    return {
        "rigid_cells": len(cells), "rigid_starts": len(poses), "ranked_poses": len(ranks),
        "baseline_outputs": len(base_rows), "baseline_poses": len(base_modes),
        "rmsds_verified": len(errors), "maximum_rmsd_discrepancy_A": max(errors),
        "single_start_rigid_fit_lower_bound_A": shapes,
        "cross_backbone_alignment": [
            {"receiver": r, "query": q, "CA_RMSD_A": ca_rmsd, "shared_ca": n}
            for (r, q), (_, ca_rmsd, n) in sorted(refs.items()) if q != r],
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("inputs", "baseline", "rigid"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    print(json.dumps(check(a.inputs, a.baseline, a.rigid), indent=2))
