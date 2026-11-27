"""Proper-fit geometry of the prepared input and the shortest ring crossing the broken bond.

This ligand-only fit measures a conformational prior. Docking success is always
measured without a ligand fit and is reported elsewhere.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from rdkit import Chem


def read_coords(path):
    return {
        int(line[6:11]): (line.split()[-1], [float(line[j:j+8]) for j in (30, 38, 46)])
        for line in path.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))
    }


def fitted_rmsd(a, b):
    x, y = a - np.mean(a, axis=0), b - np.mean(b, axis=0)
    u, _, vt = np.linalg.svd(x.T @ y)
    u[:, -1] *= np.sign(np.linalg.det(u @ vt))
    return np.mean(np.sum((x @ u @ vt - y)**2, axis=1))**0.5


def bace(inputs):
    pdbqt = inputs / "tutorial/solution/BACE_1_ligand.pdbqt"
    atoms = read_coords(pdbqt)
    mol = Chem.MolFromMol2File(str(inputs / "tutorial/data/BACE_1_ligand.mol2"), removeHs=True)
    serial = []
    for xyz in mol.GetConformer().GetPositions():
        matches = [i for i, (t, p) in atoms.items() if t not in ("HD", "G0") and np.linalg.norm(xyz - p) < 0.002]
        if len(matches) != 1:
            raise ValueError("Cannot map original MOL2 to PDBQT")
        serial.append(matches[0])
    ref = {
        line[12:16].strip(): [float(line[j:j+8]) for j in (30, 38, 46)]
        for line in (inputs / "bace/5QCT.pdb").read_text().splitlines()
        if line.startswith("HETATM") and line[17:20] == "E4Y" and line[21:26] == "A 401"
    }
    xyz = np.array([ref[a.GetProp("_TriposAtomName")] for a in mol.GetAtoms()])
    maps = mol.GetSubstructMatches(mol, uniquify=False,useChirality=True,maxMatches=10001)
    if not maps or len(maps) > 10000:
        raise ValueError("Automorphism enumeration failed")
    return ("5QCT", "crystal-derived tutorial", mol, atoms, serial, xyz[np.array(maps)])


def read_prepared(folder):
    data = json.loads((folder / "preparation.json").read_text())
    mol = Chem.MolFromSmiles(data["canonical_smiles"])
    return (data["pdb"], "CCD/ETKDG/MMFF", mol, read_coords(folder / "ligand.pdbqt"),
            data["serial"], np.array(data["mapped"]))


def measure(case):
    pdb, source, mol, atoms, serial, refs = case
    bonds = [i for i, (t, _) in atoms.items() if t == "CG0"]
    if len(bonds) != 2 or len(serial) != mol.GetNumAtoms():
        raise ValueError(f"Expected one opened bond and a complete heavy atom map for {pdb}")
    a, b = [serial.index(i) for i in bonds]
    if mol.GetBondBetweenAtoms(a, b) is None:
        raise ValueError(f"Closure carbons are not bonded in the complete graph of {pdb}")
    cut = Chem.RWMol(mol)
    cut.RemoveBond(a, b)
    ring = list(Chem.GetShortestPath(cut, a, b))
    if len(ring) < 7:
        raise ValueError(f"Unrecognized macrocycle in {pdb}")
    query = np.array([atoms[i][1] for i in serial])
    heavy_rmsd = min(fitted_rmsd(query, ref) for ref in refs)
    ring_rmsd = min(fitted_rmsd(query[ring], ref[ring]) for ref in refs)
    bound = (len(ring)/len(query))**0.5 * ring_rmsd
    return dict(pdb=pdb, source=source, heavy=len(query), ring_heavy=len(ring),
                ring_serials=" ".join(str(serial[i]) for i in ring), reference_maps=len(refs),
                whole_conformer_fit_A=f"{heavy_rmsd:.6f}", ring_fit_A=f"{ring_rmsd:.6f}",
                whole_rmsd_lower_bound_if_ring_frozen_A=f"{bound:.6f}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("historical_inputs", type=Path)
    p.add_argument("csv", type=Path)
    p.add_argument("extra_prepared", nargs="*", type=Path)
    args = p.parse_args()
    hist = args.historical_inputs
    cases = [bace(hist)] + [read_prepared(d) for d in (
        hist / "hsp/1YET", hist / "extra/1FKF", hist / "extra/1FKB", *args.extra_prepared)]
    rows = [measure(c) for c in cases]
    with args.csv.open("w", newline="") as f:
        out = csv.DictWriter(f, rows[0], lineterminator="\n")
        out.writeheader()
        out.writerows(rows)
