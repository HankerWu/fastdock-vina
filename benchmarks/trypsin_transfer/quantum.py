"""Frozen-coordinate local quantum transfer with matched water and pocket legs."""

import argparse
import hashlib
import itertools
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation
from tblite.interface import Calculator

from data import SEEDS, LIGANDS, write
from score_physical_pool import serial_index_pairs
from select_zinc_pose import read_models

BOHR = 0.529177210903
KCAL = 627.5094740631
# Each tuple selects the fragment composition, His57 hydrogen and continuum.
CASES = (
    ("acetate", "none", "water"),
    ("serine", "none", "water"),
    ("triad", "delta", "water"),
    ("triad", "epsilon", "water"),
    ("triad", "protonated", "water"),
    ("triad", "delta", "octanol"),
)


@dataclass
class Geometry:
    numbers: np.ndarray
    xyz: np.ndarray
    charge: int


def merge(first, second):
    return Geometry(
        np.concatenate((first.numbers, second.numbers)),
        np.concatenate((first.xyz, second.xyz)),
        first.charge + second.charge,
    )


def explicit_geometry(mol, xyz, mapped_h):
    """Preserve every mapped atom and place only missing H at an MMFF minimum."""
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, p)
    mol = Chem.Mol(mol)
    mol.RemoveAllConformers()
    mol.AddConformer(conf)
    full = Chem.AddHs(mol, addCoords=True)
    conf = full.GetConformer()
    fixed = set(range(mol.GetNumAtoms()))
    for parent, positions in mapped_h.items():
        hydrogens = [
            a.GetIdx()
            for a in full.GetAtomWithIdx(parent).GetNeighbors()
            if a.GetAtomicNum() == 1
        ]
        if len(hydrogens) != len(positions):
            raise ValueError("Polar hydrogen/parent identity changed")
        for idx, pos in zip(hydrogens, positions, strict=True):
            fixed.add(idx)
            conf.SetAtomPosition(idx, pos)
    mmff = AllChem.MMFFGetMoleculeProperties(full)
    if mmff is None:
        raise ValueError("No MMFF hydrogen parameters")
    ff = AllChem.MMFFGetMoleculeForceField(full, mmff)
    for i in fixed:
        ff.AddFixedPoint(i)
    if ff.Minimize(maxIts=1000):
        raise ValueError("Unresolved hydrogen geometry")
    coords = conf.GetPositions()
    if not np.allclose(coords[: len(xyz)], xyz, rtol=0, atol=1e-10):
        raise ValueError("Frozen heavy coordinates moved")
    return Geometry(
        np.array([a.GetAtomicNum() for a in full.GetAtoms()]),
        coords,
        Chem.GetFormalCharge(full),
    )


def position(line):
    return np.array([float(line[k : k + 8]) for k in (30, 38, 46)])


def ligand(lines):
    (graph,) = [
        line[14:]
        for line in lines
        if line.startswith("REMARK SMILES ")
        and not line.startswith("REMARK SMILES IDX")
    ]
    mol = Chem.MolFromSmiles(graph)
    atoms = {
        int(l[6:11]): position(l) for l in lines if l.startswith(("ATOM  ", "HETATM"))
    }
    mapped = sorted(serial_index_pairs(lines, "REMARK SMILES IDX"))
    if [i for i, _ in mapped] != list(range(1, mol.GetNumAtoms() + 1)):
        raise ValueError("Missing ligand heavy atom")
    polar = {}
    for i, s in serial_index_pairs(lines, "REMARK H PARENT"):
        polar.setdefault(i - 1, []).append(atoms[s])
    return explicit_geometry(mol, [atoms[s] for _, s in mapped], polar)


def binding_site(path):
    residues = {}
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and line[21] == "A":
            residues.setdefault((line[17:20], int(line[22:26])), {})[
                line[12:16].strip()
            ] = position(line)

    def fragment(residue, smiles, names, h=None):
        r = residues[residue]
        return explicit_geometry(
            Chem.MolFromSmiles(smiles),
            [r[n] for n in names],
            {i: [r[n]] for i, n in (h or {}).items()},
        )

    acid_names = ("CB", "CG", "OD1", "OD2")
    asp = fragment(("ASP", 189), "CC(=O)[O-]", acid_names)
    triad_asp = fragment(("ASP", 102), "CC(=O)[O-]", acid_names)
    ser = fragment(("SER", 195), "CO", ("CB", "OG"), {1: "HG"})
    hist_names = ("CB", "CG", "ND1", "CE1", "NE2", "CD2")
    his_delta = fragment(("HIS", 57), "Cc1[nH]cnc1", hist_names)
    his_epsilon = fragment(("HIS", 57), "Cc1nc[nH]c1", hist_names, {4: "HE2"})
    his_charged = fragment(("HIS", 57), "Cc1[nH]c[nH+]c1", hist_names, {4: "HE2"})
    return (
        asp,
        ser,
        {
            "delta": merge(his_delta, triad_asp),
            "epsilon": merge(his_epsilon, triad_asp),
            "protonated": merge(his_charged, triad_asp),
        },
    )


def rotations(parts, case, angles=range(0, 360, 60)):
    asp, ser, hist = parts
    group, histidine, _ = case
    if group == "acetate":
        return [asp]
    if not np.all(ser.numbers[:2] == [6, 8]):
        raise ValueError("Unexpected methanol donor")
    hid = int(
        np.where(
            (ser.numbers == 1) & (np.linalg.norm(ser.xyz - ser.xyz[1], axis=1) < 1.2)
        )[0].item()
    )
    axis = ser.xyz[1] - ser.xyz[0]
    axis /= np.linalg.norm(axis)
    geoms = []
    for phi in angles:
        xyz = ser.xyz.copy()
        xyz[hid] = xyz[1] + Rotation.from_rotvec(np.deg2rad(phi) * axis).apply(
            xyz[hid] - xyz[1]
        )
        pocket = merge(asp, Geometry(ser.numbers, xyz, 0))
        if group == "triad":
            pocket = merge(pocket, hist[histidine])
        geoms.append(pocket)
    return geoms


def energy(geom, solvent):
    calc = Calculator("GFN2-xTB", geom.numbers, geom.xyz / BOHR, charge=geom.charge)
    calc.set("verbosity", 0)
    calc.set("accuracy", 0.05)
    calc.set("max-iter", 500)
    if solvent:
        calc.add("alpb-solvation", solvent)
    return float(calc.singlepoint().get("energy")) * KCAL


def receptor_job(job):
    path, case = job
    parts = rotations(binding_site(path), case)
    return (path.stem, case), [(g, energy(g, case[-1])) for g in parts]


def pose_job(job):
    pool, pockets = job
    receptor, name, seed, arm = pool.stem.split("_")
    raw = pool.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    pose = read_models(raw.decode().splitlines())[1]
    lig = ligand(pose)
    if lig.charge != 1:
        raise ValueError("Expected +1 amidinium")
    ligand_vacuum = energy(lig, None)
    ligand_water = energy(lig, "water")
    data = []
    for case, rotamers in pockets.items():
        group, hist, solvent = case
        for index, (p, free) in enumerate(rotamers):
            complex_energy = energy(merge(p, lig), solvent)
            data.append(
                dict(
                    receptor=receptor,
                    ligand=name,
                    seed=seed,
                    arm=arm,
                    fragment=group,
                    histidine=hist,
                    solvent=solvent,
                    angle_deg=index * 60,
                    native=float(
                        next(
                            l for l in pose if l.startswith("REMARK VINA RESULT")
                        ).split()[3]
                    ),
                    ligand_vacuum=ligand_vacuum,
                    ligand_water=ligand_water,
                    pocket=free,
                    complex=complex_energy,
                    source_sha256=digest,
                )
            )
    return data


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("prepared", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--jobs", type=int, default=3)
    a = p.parse_args()
    with ProcessPoolExecutor(a.jobs) as executor:
        receptor_files = (
            (a.prepared / (r + ".pdbqt"), case)
            for r, case in itertools.product(("1BTY", "3PTB"), CASES)
        )
        compiled = dict(executor.map(receptor_job, receptor_files))
        pools = (
            (
                a.prepared / f"{r}_{l}_{s}_{arm}.pdbqt",
                {c: compiled[r, c] for c in CASES},
            )
            for r, l, s, arm in itertools.product(
                ("1BTY", "3PTB"), LIGANDS, SEEDS, ("native", "shell")
            )
        )
        rows = []
        for part in executor.map(pose_job, pools):
            rows.extend(part)
        write(a.output, rows)
