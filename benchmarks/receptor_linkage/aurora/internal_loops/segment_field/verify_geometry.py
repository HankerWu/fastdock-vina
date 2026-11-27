"""Certify the finite support of the rigid XS field on the actual docking grid."""

import argparse
import csv
import hashlib
import io
import itertools
import json
import subprocess
from pathlib import Path

import numpy as np

SEGMENTS = {"281-291": set(range(281, 292)), "303-305": {303, 304, 305}}
CUTOFF = 8.0
SPACING = 0.375
BOX_SIZE = 26.0


def receptor(path, reader):
    records = [
        line
        for line in path.read_text().splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    coords = np.array(
        [[float(line[i : i + 8]) for i in (30, 38, 46)] for line in records]
    )
    result = subprocess.run([reader, path], capture_output=True, text=True, check=True)
    typed = np.loadtxt(io.StringIO(result.stdout))
    np.testing.assert_array_equal(typed[:, :3], coords)
    return records, coords, typed[:, 3].astype(int), typed[:, 4].astype(bool)


def summarize_segment(records, coords, mask, lo, hi, grid):
    n = len(grid)
    cells = np.zeros((n - 1, n - 1, n - 1), dtype=bool)
    nodes = np.zeros((n, n, n), dtype=bool)
    by_atom = []
    for i in np.flatnonzero(mask):
        distance = np.linalg.norm(coords[i] - np.clip(coords[i], lo, hi))
        d2 = (grid - coords[i]) ** 2
        count = 0
        if distance <= CUTOFF:
            covered = (
                d2[:, 0, None, None] + d2[None, :, 1, None] + d2[None, None, :, 2]
                <= CUTOFF**2
            )
            nodes |= covered
            count = int(covered.sum())
        by_atom.append(
            {
                "atom": records[i][12:16].strip(),
                "residue": records[i][17:26].strip(),
                "min_to_box_A": float(distance),
                "min_to_grid_vertex_A": float(np.sqrt(d2.min(axis=0).sum())),
                "vertices_in_range": count,
            }
        )
    for a, b, c in itertools.product((0, 1), repeat=3):
        cells |= nodes[a : a + n - 1, b : b + n - 1, c : c + n - 1]
    closest = min(by_atom, key=lambda atom: atom["min_to_box_A"])
    return {
        "n_XS_atoms": len(by_atom),
        "min_to_box_A": closest["min_to_box_A"],
        "closest_atom": f"{closest['residue']} {closest['atom']}",
        "min_to_grid_vertex_A": min(a["min_to_grid_vertex_A"] for a in by_atom),
        "n_atoms_reaching_box": sum(a["min_to_box_A"] <= CUTOFF for a in by_atom),
        "n_atoms_reaching_vertices": sum(a["vertices_in_range"] > 0 for a in by_atom),
        "reachable_vertices": int(nodes.sum()),
        "possibly_affected_trilinear_cells": int(cells.sum()),
    }, by_atom


def run(work, raw, reader, out):
    out.mkdir(parents=True, exist_ok=True)
    center = np.array(json.loads((work / "prepared/2X6E.json").read_text())["center"])
    n_intervals = int(np.ceil(BOX_SIZE / SPACING))
    lo = center - (n_intervals * SPACING) / 2
    hi = lo + n_intervals * SPACING
    grid = lo + np.arange(n_intervals + 1)[:, None] * SPACING
    baseline = work / "prepared/2X6E.pdbqt"
    old, _, old_xs, _ = receptor(baseline, reader)
    arm_paths = {
        "1OL5_restrained": work / "template1OL5/trans_model/restrained/receptor.pdbqt",
        "1OL7_restrained": work / "template1OL7/trans_model/restrained/receptor.pdbqt",
        "1OL7_released": work / "template1OL7/trans_model/released/receptor.pdbqt",
        **{
            arm: raw / (arm + ".pdbqt")
            for arm in ("far_both", "short_only", "long_only", "both")
        },
    }
    summaries, atoms, fields = {}, [], {}
    for arm, path in arm_paths.items():
        rows, coords, xs, active = receptor(path, reader)
        assert rows[: len(old)] == old
        np.testing.assert_array_equal(xs[: len(old)], old_xs)
        residue = np.array([int(line[22:26]) for line in rows])
        new = np.arange(len(rows)) >= len(old)
        assert set(residue[new]) == set.union(*SEGMENTS.values())
        assert np.array_equal(
            active, [line.split()[-1] not in ("H", "HD") for line in rows]
        )
        # Check order too: same ordered contributing terms imply the same
        # floating-point additions for a fixed query and ligand atom type.
        distance = np.linalg.norm(coords - np.clip(coords, lo, hi), axis=1)
        fields[arm] = [
            (tuple(xyz), int(t))
            for xyz, t in zip(
                coords[active & (distance <= CUTOFF)], xs[active & (distance <= CUTOFF)]
            )
        ]
        group = {}
        for name, ids in SEGMENTS.items():
            select = new & active & np.isin(residue, list(ids))
            group[name], individuals = summarize_segment(
                rows, coords, select, lo, hi, grid
            )
            atoms.extend(dict(arm=arm, segment=name, **a) for a in individuals)
        summaries[arm] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "unchanged_historical_XS_types": len(old),
            "segments": group,
        }
    assert fields["far_both"] == fields["short_only"]
    assert fields["both"] == fields["long_only"] == fields["1OL7_released"]
    _, dry_coords, dry_types, dry_active = receptor(baseline, reader)
    dry_near = dry_active & (
        np.linalg.norm(dry_coords - np.clip(dry_coords, lo, hi), axis=1) <= CUTOFF
    )
    assert [
        (tuple(x), int(t)) for x, t in zip(dry_coords[dry_near], dry_types[dry_near])
    ] == fields["far_both"]
    result = {
        "cutoff_A": CUTOFF,
        "spacing_A": SPACING,
        "lo_A": lo.tolist(),
        "hi_A": hi.tolist(),
        "vertices": (n_intervals + 1) ** 3,
        "cells": n_intervals**3,
        "analysis": "Reach is geometric support before any map rounding, independent of ligand type; it bounds regions of possible differences, not their magnitude.",
        "ordered_reachable_XS_field_equalities": {
            "dry = far_both = short_only": len(fields["far_both"]),
            "1OL7_released = both = long_only": len(fields["both"]),
        },
        "receptors": summaries,
    }
    (out / "geometry.json").write_text(json.dumps(result, indent=2) + "\n")
    with (out / "atom_support.csv").open("w", newline="") as dest:
        writer = csv.DictWriter(dest, fieldnames=list(atoms[0]))
        writer.writeheader()
        writer.writerows(atoms)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "raw", "reader", "out"):
        p.add_argument(name, type=Path)
    args = p.parse_args()
    run(*(x.resolve() for x in (args.work, args.raw, args.reader, args.out)))
