"""Retrospective intraligand energy ablation on saved native Late pose pools."""

import argparse
import csv
import hashlib
import itertools
import math
import statistics
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import score_physical_pool as physical

TABLES = ROOT / "benchmarks/cutoff_switch"
ARMS = ("native", "inter_only", "gas_global", "screened_global", "screened_local")
SHARED = {
    "t4l": {"1L90", "186L"},
    "hcaii": {"6WQ4", "6WQ8"},
    "brd4": {"3MXF", "4OGI", "4OGJ"},
}


def load(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def panels(evidence):
    checksums = dict(
        line.split()[::-1]
        for line in (TABLES / "brd4/poses.sha256").read_text().splitlines()
    )
    for target, entries in (
        ("t4l", "mixed_t4l.csv"),
        ("hcaii", "mixed_hcaii.csv"),
        ("brd4", "brd4/cells.csv"),
    ):
        for cell in load(TABLES / entries):
            if cell["arm"] != "native":
                continue
            lig, rec, seed = cell["ligand"], cell["receptor"], cell["seed"]
            if target == "brd4":
                filename = f"{lig}_{rec}_{seed}_native.pdbqt"
                path = evidence / target / "runs" / filename
                checksum = checksums[filename]
            else:
                path = evidence / target / "mixed" / f"{rec}_{lig}_{seed}_native.pdbqt"
                checksum = cell["sha256"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != checksum:
                raise ValueError(f"Pose checksum mismatch: {path}")
            cell["path"] = path
            yield (target, lig), cell


def screened_properties(molecule):
    props = physical.mmff_properties(molecule)
    props.SetMMFFDielectricConstant(4.0)
    return props


def screened_reference(ligand):
    """Use the historical conformer count, embedding and minimizer with epsilon=4."""
    molecule = Chem.Mol(ligand)
    molecule.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = physical.GAS_EMBED_SEED
    params.numThreads = 1
    conformers = AllChem.EmbedMultipleConfs(
        molecule, numConfs=physical.GAS_CONFORMERS, params=params
    )
    props = screened_properties(molecule)
    energies = []
    for conf in conformers:
        try:
            energies.append(physical.minimized_energy(molecule, props, conf, ()))
        except ValueError:
            continue
    if not energies:
        raise ValueError("No screened reference conformer converged")
    return min(energies), len(energies)


def fixed_energy(molecule, props):
    conf = molecule.GetConformer()
    heavy = [a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum() > 1]
    xyz = conf.GetPositions()[heavy]
    energy = physical.minimized_energy(molecule, props, conf.GetId(), heavy)
    if not np.array_equal(xyz, conf.GetPositions()[heavy]):
        raise ValueError("Frozen heavy coordinates changed")
    return energy


def pose_energies(molecule):
    gas_mol = Chem.Mol(molecule)
    gas_pose = fixed_energy(gas_mol, physical.mmff_properties(gas_mol))
    props = screened_properties(molecule)
    screen_pose = fixed_energy(molecule, props)
    screen_basin = physical.minimized_energy(
        molecule, props, molecule.GetConformer().GetId(), ()
    )
    if screen_basin > screen_pose + 1e-5:
        raise ValueError("Releasing coordinates increased internal energy")
    return gas_pose, screen_pose, screen_basin


def ligand_rows(job):
    (target, lig), cells, binary, evidence = job
    prepared = evidence / target / "prepared"
    receptor = prepared / f"{cells[0]['receptor']}_receptor.pdbqt"
    ligand = prepared / f"{lig}_ligand.pdbqt"
    scale = physical.evaluate_coordinates(binary, ligand, receptor, None, 60)["scale"]
    first = physical.read_poses(cells[0]["path"])[0]
    reference_ligand, identity = physical.decode_ligand(first)
    gas_ref, gas_n = physical.gas_reference(reference_ligand)
    screen_ref, screen_n = screened_reference(reference_ligand)
    outputs = []
    for cell in cells:
        poses = physical.read_poses(cell["path"])
        physical.fixed_coordinate_inputs(poses, ligand, None)
        rmsds = list(map(float, cell["pool_rmsds"].split(",")))
        if (
            poses[0].native_score != float(cell["score"])
            or len({p.native_unbound for p in poses}) != 1
        ):
            raise ValueError(f"Inconsistent native score in {cell['path']}")
        for pose, rmsd in zip(poses, rmsds, strict=True):
            molecule, current = physical.decode_ligand(pose)
            if current != identity:
                raise ValueError(f"Chemical identity changed at {cell['path']}")
            gas_pose, screen_pose, screen_basin = pose_energies(molecule)
            inter_only = pose.native_score - scale * (
                pose.native_intra - pose.native_unbound
            )
            weight = scale * physical.STRAIN_WEIGHT
            outputs.append(
                {
                    "target": target,
                    "ligand": lig,
                    "receptor": cell["receptor"],
                    "seed": int(cell["seed"]),
                    "model": pose.model,
                    "measurement": float(cell.get("measurement", cell.get("kd_uM"))),
                    "rmsd": rmsd,
                    "scale": scale,
                    "native_intra": pose.native_intra,
                    "native_unbound": pose.native_unbound,
                    "eps1_fixed": gas_pose,
                    "eps1_reference": gas_ref,
                    "eps1_conformers": gas_n,
                    "eps4_fixed": screen_pose,
                    "eps4_reference": screen_ref,
                    "eps4_conformers": screen_n,
                    "eps4_basin": screen_basin,
                    "native": pose.native_score,
                    "inter_only": inter_only,
                    "gas_global": inter_only + weight * (gas_pose - gas_ref),
                    "screened_global": inter_only + weight * (screen_pose - screen_ref),
                    "screened_local": inter_only
                    + weight * (screen_pose - screen_basin),
                }
            )
    return outputs


def evaluate(rows):
    pools = defaultdict(list)
    for r in rows:
        pools[(r["target"], r["receptor"], r["seed"], r["ligand"])].append(r)
    cells, pairs = [], []
    for (target, receptor, seed, ligand), pool in sorted(pools.items()):
        for arm in ARMS:
            ranked = sorted(pool, key=lambda row: (row[arm], row["model"]))
            winner = ranked[0]
            near = next((i for i, r in enumerate(ranked, 1) if r["rmsd"] <= 2), None)
            cells.append(
                {
                    "target": target,
                    "receptor": receptor,
                    "seed": seed,
                    "ligand": ligand,
                    "arm": arm,
                    "measurement": winner["measurement"],
                    "model": winner["model"],
                    "score": winner[arm],
                    "rmsd": winner["rmsd"],
                    "first_near_native": near,
                    "best_rmsd": min(r["rmsd"] for r in pool),
                }
            )
    blocks = defaultdict(list)
    for cell in cells:
        if cell["receptor"] in SHARED[cell["target"]]:
            blocks[
                (cell["target"], cell["receptor"], cell["seed"], cell["arm"])
            ].append(cell)
    for (target, receptor, seed, arm), cells_in_block in sorted(blocks.items()):
        for strong, weak in itertools.combinations(
            sorted(cells_in_block, key=lambda r: r["measurement"]), 2
        ):
            ratio = weak["measurement"] / strong["measurement"]
            if ratio >= 3:
                gap = weak["score"] - strong["score"]
                pairs.append(
                    {
                        "target": target,
                        "receptor": receptor,
                        "seed": seed,
                        "arm": arm,
                        "strong": strong["ligand"],
                        "weak": weak["ligand"],
                        "measurement_ratio": ratio,
                        "weak_minus_strong": gap,
                        "credit": float(gap > 0) + 0.5 * (gap == 0),
                    }
                )
    summary = []
    for target, arm in itertools.product(SHARED, ARMS):
        cell = [r for r in cells if r["target"] == target and r["arm"] == arm]
        pair = [r for r in pairs if r["target"] == target and r["arm"] == arm]
        summary.append(
            {
                "target": target,
                "arm": arm,
                "cells": len(cell),
                "top1": sum(r["rmsd"] <= 2 for r in cell),
                "median_rmsd": statistics.median(r["rmsd"] for r in cell),
                "near_in_pool": sum(r["first_near_native"] is not None for r in cell),
                "med_rank_if_available": statistics.median(
                    [r["first_near_native"] for r in cell if r["first_near_native"]]
                    or [math.nan]
                ),
                "pairs": len(pair),
                "credit": sum(r["credit"] for r in pair),
                "unique_pairs": len({(r["strong"], r["weak"]) for r in pair}),
            }
        )
    return cells, pairs, summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "evidence", type=Path, help="Extracted FastDock-Vina-Pro-results.zip root"
    )
    p.add_argument("binary", type=Path, help="Vina with net-polar normalization output")
    p.add_argument("output", type=Path)
    p.add_argument("--workers", type=int, default=2)
    args = p.parse_args()
    by_ligand = defaultdict(list)
    for key, cell in panels(args.evidence):
        by_ligand[key].append(cell)
    jobs = [
        (key, cells, args.binary.resolve(), args.evidence)
        for key, cells in sorted(by_ligand.items())
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = [r for group in pool.map(ligand_rows, jobs) for r in group]
    cells, pairs, summary = evaluate(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    for name, data in (
        ("poses", rows),
        ("cells", cells),
        ("pairs", pairs),
        ("summary", summary),
    ):
        write(args.output / (name + ".csv"), data)


if __name__ == "__main__":
    main()
