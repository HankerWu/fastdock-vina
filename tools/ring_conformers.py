"""Prepare a matched pair of ligand PDBQT inputs with different small-ring geometries."""

import argparse
import json
from pathlib import Path

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem


def fit(mobile, fixed):
    """Return a proper rotation and translation in the row-vector convention."""
    u, _, v = np.linalg.svd((mobile - mobile.mean(0)).T @ (fixed - fixed.mean(0)))
    rot = u @ np.diag([1.0, 1.0, np.linalg.det(u @ v)]) @ v
    return rot, fixed.mean(0) - mobile.mean(0) @ rot


def ring_fragments(mol):
    fragments = []
    for ring in mol.GetRingInfo().AtomRings():
        if not 5 <= len(ring) <= 6 or all(
            mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring
        ):
            continue
        fragment = set(ring)
        for i in ring:
            fragment.update(
                a.GetIdx()
                for a in mol.GetAtomWithIdx(i).GetNeighbors()
                if a.GetAtomicNum() > 1
            )
        fragments.append(sorted(fragment))
    return fragments


def ring_distance(first, second, fragments, permutations):
    """Largest ring/first-shell proper-fit RMSD, minimized over graph symmetries."""
    if not fragments:
        return 0.0
    best = np.inf
    for perm in permutations:
        errors = []
        for fragment in fragments:
            a = first[fragment]
            b = second[np.asarray(perm)[fragment]]
            rot, shift = fit(a, b)
            errors.append(np.sqrt(np.mean(np.sum((a @ rot + shift - b) ** 2, axis=1))))
        best = min(best, max(errors))
    return float(best)


def prepare_pair(mol, folder, *, seed=20261126, count=64, window=12.0, distance=0.35):
    """Ignore supplied coordinates and retain stereochemistry and ionization of the input graph."""
    folder.mkdir(parents=True, exist_ok=True)
    mol = Chem.RemoveHs(mol)
    if len(Chem.GetMolFrags(mol)) != 1 or any(
        a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        raise ValueError("Expected one connected, non-radical small molecule")
    mol.RemoveAllConformers()
    fragments = ring_fragments(mol)
    perms = mol.GetSubstructMatches(
        mol, uniquify=False, useChirality=True, maxMatches=1025
    )
    if len(perms) > 1024:
        raise ValueError("More than 1024 ligand automorphisms")
    full = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1
    params.useSmallRingTorsions = True
    ids = list(AllChem.EmbedMultipleConfs(full, numConfs=count, params=params))
    if not ids or not AllChem.MMFFHasAllMoleculeParams(full):
        raise ValueError("ETKDG embedding failed or MMFF parameters unavailable")
    results = AllChem.MMFFOptimizeMoleculeConfs(
        full, numThreads=1, mmffVariant="MMFF94s", maxIters=1200
    )
    ranked = sorted(
        (e, i)
        for i, (status, e) in zip(ids, results, strict=True)
        if status == 0 and np.isfinite(e)
    )
    if not ranked:
        raise ValueError("No converged conformer")
    e0, i0 = ranked[0]
    e1, i1, delta = e0, i0, 0.0
    xyz = full.GetConformer(i0).GetPositions()
    for e, i in ranked[1:]:
        if e - e0 > window or not fragments:
            break
        d = ring_distance(xyz, full.GetConformer(i).GetPositions(), fragments, perms)
        if d >= distance:
            e1, i1, delta = e, i, d
            break
    # Negative control: change the full input conformer while retaining its ring neighborhoods.
    sham_id, sham_energy = i0, e0
    for e, i in ranked[1:]:
        if e - e0 > window:
            break
        candidate = full.GetConformer(i).GetPositions()
        if ring_distance(xyz, candidate, fragments, perms) > 0.15:
            continue
        if (
            ring_distance(xyz, candidate, [list(range(mol.GetNumAtoms()))], perms)
            >= 0.5
        ):
            sham_id, sham_energy = i, e
            break
    first_e, first_id = min(ranked, key=lambda entry: entry[1])
    topologies = []
    for label, cid in [
        ("base", i0),
        ("alternate", i1),
        ("sham", sham_id),
        ("first", first_id),
    ]:
        conformer = Chem.Mol(full)
        conformer.RemoveAllConformers()
        conformer.AddConformer(full.GetConformer(cid), assignId=True)
        (setup,) = MoleculePreparation().prepare(conformer)
        qt, ok, error = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
        if not ok:
            raise ValueError(error)
        topologies.append(
            [
                line[:30] + line[54:] if line.startswith(("ATOM", "HETATM")) else line
                for line in qt.splitlines()
            ]
        )
        (folder / f"{label}.pdbqt").write_text(qt)
    if any(t != topologies[0] for t in topologies[1:]):
        raise ValueError(
            "Conformers produce unequal PDBQT torsion trees, charges or atom types"
        )
    info = {
        "seed": seed,
        "embedded": len(ids),
        "converged": len(ranked),
        "ids": [i0, i1],
        "energies_MMFF94s": [e0, e1],
        "ring_distance_A": delta,
        "fragments": fragments,
        "min_distance_A": distance,
        "energy_window": window,
        "smiles": Chem.MolToSmiles(mol),
        "alternate_distinct": i0 != i1,
        "sham_distinct": i0 != sham_id,
        "sham_energy_MMFF94s": sham_energy,
        "sham_id": sham_id,
        "first_id": first_id,
        "first_energy_MMFF94s": first_e,
        "first_distinct": first_id != i0,
    }
    (folder / "states.json").write_text(json.dumps(info, indent=2) + "\n")
    return info


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smiles",
        required=True,
        help="One fully specified ligand, including charge and stereochemistry",
    )
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    parsed = Chem.MolFromSmiles(args.smiles)
    if parsed is None:
        parser.error("Invalid SMILES")
    print(json.dumps(prepare_pair(parsed, args.out), indent=2))
