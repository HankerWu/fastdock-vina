"""Native ranked modes in a single receptor frame, with alternative site occupancy."""

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from scipy.spatial import cKDTree

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

from rank_halogen_pool import models, pdb_coordinates, rank_pool
from run_assembly import PROBES, CONTROLS, INITIAL, REPEAT


def independent_positions(folder, destination, ref):
    """Verify transformed receptor coordinates and identify matching official ligand copies by atom names."""
    cc = folder.name.split("_")[1]
    original = [
        l
        for l in (folder / "deposition.pdb").read_text().splitlines()
        if l.startswith("HETATM")
        and l[17:20] == cc
        and l[16] in " A"
        and l[76:78].strip() not in ("H", "D")
    ]
    rx = ref.GetConformer().GetPositions()
    nearest, lookup = cKDTree(np.array([pdb_coordinates(l) for l in original])).query(
        rx
    )
    if nearest.max() > 0.01:
        raise ValueError("mismatched PDB and SDF coordinates")
    names = [original[i][12:16] for i in lookup]
    with gzip.open(destination / "deposited-assembly.pdb.gz", "rt") as f:
        deposit = f.read().splitlines()
    models = []
    group = []
    for l in deposit:
        if l.startswith("MODEL") and group:
            models.append(group)
            group = []
        group.append(l)
        if l.startswith("ENDMDL"):
            models.append(group)
            group = []
    if group:
        models.append(group)
    copies = []
    for model in models:
        byres = {}
        for l in model:
            if l.startswith("HETATM") and l[17:20] == cc and l[16] in " A":
                byres.setdefault(l[21:27], {})[l[12:16]] = l
        for atoms in byres.values():
            if not all(k in atoms for k in names):
                continue
            xyz = np.array([pdb_coordinates(atoms[n]) for n in names])
            if np.linalg.norm(xyz.mean(0) - rx.mean(0)) <= 2:
                copies.append(xyz)
    if not copies:
        raise ValueError(f"no colocated official reference for {cc}")
    # All included receptor heavy atoms must have a matching deposited atom of the same identity and position.
    dep_groups = {}
    for l in deposit:
        if (
            l.startswith("ATOM")
            and l[16] in " A"
            and l[76:78].strip() not in ("H", "D")
        ):
            dep_groups.setdefault(
                (l[12:16].strip(), l[17:20], int(l[22:26])), []
            ).append(pdb_coordinates(l))
    dep_trees = {k: cKDTree(v) for k, v in dep_groups.items()}
    terminal_keys = {(a, b) for n, a, b in dep_groups if n == "OXT"}
    for a, b in terminal_keys:
        if ("O", a, b) in dep_groups:
            shared = cKDTree(dep_groups[("O", a, b)] + dep_groups[("OXT", a, b)])
            dep_trees[("O", a, b)] = shared
            dep_trees[("OXT", a, b)] = shared
    max_error = 0.0
    total = 0
    included = []
    for l in (destination / "receptor.pdbqt").read_text().splitlines():
        if l.startswith("ATOM") and l.split()[-1] not in ("HD", "H"):
            key = (l[12:16].strip(), l[17:20], int(l[22:26]))
            xyz = pdb_coordinates(l)
            if key not in dep_trees:
                raise ValueError(f"unexpected heavy atom {key}")
            max_error = max(max_error, float(dep_trees[key].query(xyz)[0]))
            total += 1
            included.append(xyz)
    if max_error > 0.005:
        raise ValueError(f"assembly coordinate error {max_error}")
    inside = np.vstack([v for v in dep_groups.values()])
    represented = cKDTree(included).query(inside)[0] < 0.005
    missing = inside[~represented]
    nearest_missing = (
        float(cKDTree(np.vstack(copies)).query(missing)[0].min())
        if len(missing)
        else None
    )
    if nearest_missing is not None and nearest_missing < 6:
        raise ValueError(f"near-site missing assembly atom {nearest_missing}")
    return copies, dict(
        verified_atoms=total,
        max_official_coordinate_error_A=max_error,
        same_site_copies=len(copies),
        closest_missing_atom_A=nearest_missing,
    )


def run(work, dest, out, high=False):
    rows = []
    cases = []
    modes = []
    for case in (PROBES + ("2ROX_T44",)) if high else PROBES + CONTROLS:
        folder = work / case
        directory = dest / case
        ref = Chem.SDMolSupplier(str(folder / "crystal.sdf"), removeHs=True)[0]
        copies, info = independent_positions(folder, directory, ref)
        refs = [ref]
        for xyz in copies:
            dup = Chem.Mol(ref)
            dup.GetConformer().SetPositions(xyz)
            refs.append(dup)
        receptor_lines = (directory / "receptor.pdbqt").read_text().splitlines()
        added_chains = {
            c["chain"]
            for c in json.loads((directory / "assembly.json").read_text())["copies"]
            if c["chain"] != c["source"]
        }
        added = [
            pdb_coordinates(l)
            for l in receptor_lines
            if l.startswith("ATOM")
            and l[21] in added_chains
            and l.split()[-1] not in ("HD", "H")
        ]
        number = (
            int((cKDTree(added).query(ref.GetConformer().GetPositions())[0] <= 4).sum())
            if added
            else 0
        )
        cases.append(
            dict(case=case, added_chain_contacts_ligand_atoms_4A=number, **info)
        )
        for seed in (
            (316227, 331663, 346410)
            if high
            else INITIAL + (REPEAT if case in PROBES else ())
        ):
            for arm in ("asu", "assembly"):
                for cap in (60000, 240000) if high else (60000,):
                    path = (
                        directory / f"{arm}-{seed}-cap{cap}.pdbqt"
                        if high
                        else folder / f"poses-{seed}.pdbqt"
                        if arm == "asu" and seed in INITIAL
                        else directory / f"{arm}-{seed}.pdbqt"
                    )
                    distances = []
                    rr = (folder if arm == "asu" else directory) / "receptor.pdbqt"
                    placement = rank_pool(rr, path)
                    for idx, score, mol, xyz in models(path):
                        conf = Chem.Conformer(mol.GetNumAtoms())
                        conf.SetPositions(xyz)
                        mol.AddConformer(conf)
                        rms = [
                            rdMolAlign.CalcRMS(mol, r, maxMatches=100000) for r in refs
                        ]
                        distances.append((rms[0], min(rms)))
                        modes.append(
                            dict(
                                case=case,
                                seed=seed,
                                arm=arm,
                                cap=cap,
                                model=idx,
                                score=score,
                                rmsd_A=rms[0],
                                site_rmsd_A=min(rms),
                            )
                        )
                    rows.append(
                        dict(
                            case=case,
                            seed=seed,
                            arm=arm,
                            cap=cap,
                            top1_A=distances[0][0],
                            top1_site_A=distances[0][1],
                            oracle_A=min(v[0] for v in distances),
                            oracle_site_A=min(v[1] for v in distances),
                            directional_top1_site_A=distances[
                                next(
                                    i
                                    for i, v in enumerate(placement)
                                    if v["directional_rank"] == 1
                                )
                            ][1],
                            radial_top1_site_A=distances[
                                next(
                                    i
                                    for i, v in enumerate(placement)
                                    if v["radial_rank"] == 1
                                )
                            ][1],
                        )
                    )
    for name, data in [
        ("assembly-cells", rows),
        ("assembly-modes", modes),
        ("assembly-audit", cases),
    ]:
        with (out / (name + ("-high" if high else "") + ".csv")).open("w") as f:
            w = csv.DictWriter(f, fieldnames=data[0])
            w.writeheader()
            w.writerows(data)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("dest", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--high", action="store_true")
    a = p.parse_args()
    run(a.work, a.dest, a.out, a.high)
