"""Independent pose correspondence, geometry checks, and matched candidate-stage audit."""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from prepare_holdout import PRIMARY, ALTERNATES, POSE

sys.path.insert(0, str(POSE))
import audit
import paired


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def frames(work, cases):
    data = {pdb: audit.reference(work, pdb, *cases[pdb][1:]) for pdb in cases}
    targets, inputs, lower = {}, {}, {}
    for q in cases:
        mol, xyz, tags = audit.coordinates(
            (work / q / "ligand.pdbqt").read_text().splitlines()
        )
        pairs = np.array(
            [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
        )
        lengths = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1)
        inputs[q] = (tags, pairs, lengths)
        qc, center, qmol = data[q]
        ref = qmol.GetConformer().GetPositions()
        matches = qmol.GetSubstructMatches(
            mol, uniquify=False, useChirality=True, maxMatches=10001
        )
        if not 0 < len(matches) < 10001:
            raise ValueError(q)
        lower[q] = min(
            Rotation.align_vectors(ref[list(m)] - ref.mean(0), xyz - xyz.mean(0))[1]
            / np.sqrt(len(xyz))
            for m in matches
        )
        for r in cases:
            rc, rc_center, _ = data[r]
            keys = sorted(
                k
                for k in qc.keys() & rc.keys()
                if np.linalg.norm(qc[k] - center) < 14
                and np.linalg.norm(rc[k] - rc_center) < 14
            )
            x, y = np.array([qc[k] for k in keys]), np.array([rc[k] for k in keys])
            rot = Rotation.align_vectors(y - y.mean(0), x - x.mean(0))[0]
            targets[q, r] = qmol, rot.apply(ref - x.mean(0)) + y.mean(0)
    return targets, inputs, lower


def main(work):
    cases = {
        p: (PRIMARY | ALTERNATES)[p]
        for p in json.loads((work / "cases.json").read_text())
    }
    targets, inputs, lower = frames(work, cases)
    expected = {
        (r["key"], int(r["mode"])): float(r["rmsd_A"]) for r in read(work / "modes.csv")
    }
    largest, bond_error, inversions, checked = 0.0, 0.0, 0, 0
    geom = []
    for key in sorted({k for k, i in expected}):
        q, r, *_ = key.split("_")
        tags, pairs, lengths = inputs[q]
        crystal, reference = targets[q, r]
        for i, block in enumerate(
            (work / "runs" / (key + ".pdbqt")).read_text().split("ENDMDL")[:-1], 1
        ):
            mol, xyz, detected = audit.coordinates(block.splitlines())
            flips = sum(
                t != "CHI_UNSPECIFIED" and detected[j] != t for j, t in tags.items()
            )
            diff = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1) - lengths
            matches = crystal.GetSubstructMatches(
                mol, uniquify=False, useChirality=True, maxMatches=10001
            )
            if not 0 < len(matches) < 10001:
                raise ValueError((key, i))
            rms = min(
                np.linalg.norm(reference[list(m)] - xyz) / np.sqrt(len(xyz))
                for m in matches
            )
            error = abs(rms - expected.pop((key, i)))
            geom.append(
                dict(
                    key=key,
                    mode=i,
                    independent_rmsd_A=rms,
                    disagreement_A=error,
                    max_bond_change_A=max(abs(diff)),
                    inverted_centers=flips,
                )
            )
            largest = max(largest, error)
            bond_error = max(bond_error, max(abs(diff)))
            inversions += bool(flips)
            checked += 1
    if expected or largest > 1e-5 or inversions or bond_error > 0.005:
        raise ValueError((len(expected), largest, inversions, bond_error))
    summaries = {}
    stages = read(work / "stages.csv")
    all_rows = {(r["key"], r["stage"], r["input_index"]): r for r in stages}
    for stage in ("search", "before", "after"):
        paired_rows = [
            (
                r,
                all_rows[
                    r["key"].removesuffix("linear") + "shift", stage, r["input_index"]
                ],
            )
            for r in stages
            if r["key"].endswith("_linear") and r["stage"] == stage
        ]
        summaries[stage] = dict(
            comparisons=len(paired_rows),
            identical_coords=sum(
                a["coords_sha256"] == b["coords_sha256"] for a, b in paired_rows
            ),
            identical_coords_and_energy=sum(
                a["sha256"] == b["sha256"] for a, b in paired_rows
            ),
            linear_hits=sum(float(a["rmsd_A"]) <= 2 for a, b in paired_rows),
            shift_hits=sum(float(b["rmsd_A"]) <= 2 for a, b in paired_rows),
        )
    if not all(
        summaries[k]["comparisons"] == summaries[k]["identical_coords_and_energy"]
        for k in ("before", "search")
    ):
        raise ValueError("Grid-stage outputs differ")
    for arm in ("linear", "shift"):
        changes = defaultdict(int)
        for r in stages:
            if r["stage"] != "before" or not r["key"].endswith("_" + arm):
                continue
            after = all_rows[r["key"], "after", r["input_index"]]
            before_hit, after_hit = float(r["rmsd_A"]) <= 2, float(after["rmsd_A"]) <= 2
            changes[f"{int(before_hit)}_to_{int(after_hit)}"] += 1
        summaries[arm + "_events"] = changes
    paired.write(work / "geometry.csv", geom)
    cells = read(work / "cells.csv")
    by_key = {row["key"]: row for row in cells}
    if len(by_key) != len(cells):
        raise ValueError("Duplicate docking case")
    counters = (
        "hunt_calls",
        "hunt_evals",
        "second_calls",
        "second_evals",
        "explicit_inputs",
    )
    for row in cells:
        key = row["key"]
        if row["arm"] not in ("linear", "shift"):
            continue
        before = [r for r in stages if r["key"] == key and r["stage"] == "before"]
        after = [r for r in stages if r["key"] == key and r["stage"] == "after"]
        if len(before) != len(after) or len(before) != int(row["explicit_inputs"]):
            raise ValueError("Stage trace count mismatch for " + key)
        baseline = by_key[key.rsplit("_", 1)[0] + "_linear"]
        if any(row[field] != baseline[field] for field in counters):
            raise ValueError("Grid-stage counter mismatch for " + key)
    # Matched rigid-receptor compatibility of the experimental poses, including
    # the overlay assessment for a ligand without a crystal in that receiver.
    reference_packing = {}
    for (query, receiver), (crystal_mol, xyz) in targets.items():
        raw = (work / receiver / "receptor.pdbqt").read_text().splitlines()
        rec = list(paired.atoms(raw).values())
        dist = np.linalg.norm(
            xyz[:, None, :] - np.array([a[0] for a in rec])[None, :, :], axis=2
        )
        rad = np.array([paired.RADII[a.GetSymbol()] for a in crystal_mol.GetAtoms()])
        cutoff = 0.75 * (rad[:, None] + np.array([a[1] for a in rec])[None, :])
        reference_packing[query + "_" + receiver] = {
            "clash_pairs": int((dist < cutoff).sum()),
            "minimum_heavy_distance_A": float(dist.min()),
        }
    data = dict(
        reference_packing=reference_packing,
        modes_checked=checked,
        max_rmsd_disagreement_A=largest,
        max_bond_change_A=bond_error,
        poses_with_chiral_changes=inversions,
        prepared_to_aligned_crystal_rmsd_A=lower,
        stages=summaries,
    )
    (work / "verification.json").write_text(json.dumps(data, indent=2) + "\n")
    print(json.dumps(data, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    main(p.parse_args().work)
