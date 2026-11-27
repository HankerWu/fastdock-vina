"""Independent protein-frame RMSD, coordinate chemistry and paired CLI checks."""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import check_predictions as previous
from frame_check import positions


def serial_map(path):
    mapping = {}
    for row in path.read_text().splitlines():
        if row.startswith("REMARK INDEX MAP "):
            values = [int(x) for x in row.split()[3:]]
            mapping.update(zip(values[::2], values[1::2], strict=True))
    return mapping


def counters(path):
    markers = ("Search trajectory:", "Hunt chain:", "Pocket start:", "Search counters:")
    return [
        line
        for line in path.with_suffix(".log").read_text().splitlines()
        if line.startswith(markers)
    ]


def chemistry(work, drug):
    path = work / "prepared" / f"{drug}.pdbqt"
    mol = Chem.MolFromMolFile(str(path.with_suffix(".sdf")), removeHs=False)
    mapping = serial_map(path)
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    source = mol.GetConformer().GetPositions()[heavy]
    serials = [mapping[i + 1] for i in heavy]
    inverse = {i: j for j, i in enumerate(heavy)}
    bonds = np.array(
        [
            [inverse[b.GetBeginAtomIdx()], inverse[b.GetEndAtomIdx()]]
            for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in inverse and b.GetEndAtomIdx() in inverse
        ]
    )
    # Signed volume uses a fixed order of three heavy neighbours; the implicit
    # fourth substituent is H for the specified stereocentre in this data set.
    centres = []
    for i, tag in Chem.FindMolChiralCenters(mol, includeUnassigned=True):
        if tag == "?":
            raise ValueError(f"Unspecified stereocentre {drug} {i}")
        heavy_neighbors = [
            inverse[a.GetIdx()]
            for a in mol.GetAtomWithIdx(i).GetNeighbors()
            if a.GetIdx() in inverse
        ]
        if len(heavy_neighbors) != 3:
            raise ValueError("This check expects a three-heavy-neighbour stereocentre")
        centres.append((inverse[i], heavy_neighbors))
    rings = [
        tuple(inverse[i] for i in ring)
        for ring in mol.GetRingInfo().AtomRings()
        if len(ring) == 3
    ]
    return source, serials, bonds, centres, rings


def lengths(x, pairs):
    return np.linalg.norm(x[pairs[:, 1]] - x[pairs[:, 0]], axis=1)


def run(work, runs):
    independent_rmsd = previous.audit_blind(work, [runs / "historical", runs / "fresh"])
    chemical = {}
    for drug in ("at9283", "danusertib"):
        ref, serials, bonds, centres, rings = chemistry(work, drug)
        ring_edges = np.array(
            [(r[k], r[(k + 1) % 3]) for r in rings for k in range(3)], dtype=int
        )
        reference_lengths = lengths(ref, bonds)
        expected_sign = np.sign(
            [np.linalg.det(ref[neighbors] - ref[i]) for i, neighbors in centres]
        )
        max_bond = 0.0
        max_ring = 0.0
        ring_max = -np.inf
        ring_min = np.inf
        inverted = 0
        modes = 0
        for file in sorted(runs.glob(f"*/blind/*_{drug}_*.pdbqt")):
            for pose in positions(file):
                xyz = np.array([pose[s] for s in serials])
                if not np.all(np.isfinite(xyz)):
                    raise ValueError(f"Nonfinite pose: {file}")
                max_bond = max(
                    max_bond,
                    float(np.max(np.abs(lengths(xyz, bonds) - reference_lengths))),
                )
                sign = np.sign(
                    [np.linalg.det(xyz[neighbors] - xyz[i]) for i, neighbors in centres]
                )
                inverted += int(np.count_nonzero(sign != expected_sign))
                if len(ring_edges):
                    d = lengths(xyz, ring_edges)
                    max_ring = max(
                        max_ring, float(np.max(np.abs(d - lengths(ref, ring_edges))))
                    )
                    ring_min = min(ring_min, float(d.min()))
                    ring_max = max(ring_max, float(d.max()))
                modes += 1
        # SDF at four decimals, PDBQT input/output at three; rigid transformations
        # preserve true distances. This conservative bound includes both roundings.
        assert max_bond < 0.004 and inverted == 0
        chemical[drug] = {
            "finite_modes": modes,
            "sdf_heavy_bonds_per_mode": len(bonds),
            "max_bond_length_change_vs_SDF_A": max_bond,
            "specified_stereocentres_per_mode": len(centres),
            "inverted_centres": inverted,
            "three_membered_rings_per_mode": len(rings),
            "ring_edge_range_A": [ring_min, ring_max] if rings else None,
            "max_triangle_edge_change_vs_SDF_A": max_ring if rings else None,
        }
    pair_counts = defaultdict(int)
    eval_checks = 0
    for group in ("historical", "fresh"):
        result = runs / group
        for row in previous.rows(result / "cells.csv"):
            arm, drug, seed = row["arm"], row["drug"], row["seed"]
            f = result / "blind" / f"{arm}_{drug}_{seed}.pdbqt"
            log = counters(f)
            if arm != "official":
                stats = dict()
                for line in log:
                    for term, count in re.findall(
                        r"\b(hunt_evals|second_evals|explicit_evals)=(\d+)", line
                    ):
                        stats[term] = stats.get(term, 0) + int(count)
                assert set(stats) == {"hunt_evals", "second_evals", "explicit_evals"}
                assert sum(stats.values()) == int(row["evals"])
                eval_checks += 1
            expected = []
            if arm == "far_both":
                expected.append(
                    (
                        "short vs neither",
                        result / "blind" / f"short_only_{drug}_{seed}.pdbqt",
                    )
                )
            if arm == "both":
                expected.append(
                    (
                        "both vs long",
                        result / "blind" / f"long_only_{drug}_{seed}.pdbqt",
                    )
                )
            if group == "historical" and arm in ("far_both", "both"):
                label = "dry" if arm == "far_both" else "1OL7_released"
                source = work / ("base" if arm == "far_both" else "templates")
                expected.append(
                    (
                        "historical source",
                        source / "blind" / f"{label}_{drug}_{seed}.pdbqt",
                    )
                )
            if arm.endswith("_graph") and drug == "danusertib":
                expected.append(
                    (
                        "danusertib native vs explicit",
                        result / "blind" / f"{arm[:-6]}_{drug}_{seed}.pdbqt",
                    )
                )
            for key, other in expected:
                assert f.read_bytes() == other.read_bytes()
                assert log and log == counters(other)
                pair_counts[key] += 1
    print(
        json.dumps(
            {
                "protein_frame_RMSD": independent_rmsd,
                "coordinate_chemistry": chemical,
                "CLI_counter_checks": eval_checks,
                "identical_pool_and_counter_pairs": dict(pair_counts),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("runs", type=Path)
    args = p.parse_args()
    run(args.work.resolve(), args.runs.resolve())
