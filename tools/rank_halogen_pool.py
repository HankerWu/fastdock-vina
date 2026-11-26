"""Rerank a single-ligand fixed Meeko/Vina pose pool for carbonyl-directed halogen bonds."""

import argparse
import csv
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem

# vdW sums from Auffinger et al. (PNAS 2004); weights from Koebel et al. (2016).
PARAMETERS = {"Cl": (3.27, 0.265), "Br": (3.37, 0.32), "I": (3.50, 0.40)}
WIDTH = 0.35
AMINO_ACIDS = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
    "ASH GLH HID HIE HIP CYM CYX LYN MSE SEC".split()
)


def pdb_coordinates(line):
    return [float(line[i : i + 8]) for i in (30, 38, 46)]


def carbonyls(receptor):
    coordinates = []
    lines = receptor.read_text().splitlines()
    terminal = {
        line[17:27]
        for line in lines
        if line.startswith(("ATOM  ", "HETATM")) and line[12:16].strip() == "OXT"
    }
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")) or line.split()[-1] != "OA":
            continue
        atom = line[12:16].strip()
        residue = line[17:20].strip()
        if (residue in AMINO_ACIDS and atom == "O" and line[17:27] not in terminal) or (
            residue,
            atom,
        ) in {
            ("ASN", "OD1"),
            ("GLN", "OE1"),
        }:
            coordinates.append(pdb_coordinates(line))
    return np.array(coordinates, dtype=float).reshape((-1, 3))


def donors(mol):
    for atom in mol.GetAtoms():
        element = atom.GetSymbol()
        neighbors = list(atom.GetNeighbors())
        if (
            element in PARAMETERS
            and len(neighbors) == 1
            and neighbors[0].GetAtomicNum() == 6
        ):
            yield atom.GetIdx(), neighbors[0].GetIdx(), element


def corrections(mol, coordinates, acceptors):
    """An acceptor can compete for each donor; one oxygen contributes per halogen."""
    distance_only = directional = 0.0
    if len(acceptors) == 0:
        return distance_only, directional
    for halogen, carbon, element in donors(mol):
        center, weight = PARAMETERS[element]
        axis = coordinates[carbon] - coordinates[halogen]
        to_oxygen = acceptors - coordinates[halogen]
        length = np.linalg.norm(to_oxygen, axis=1)
        if np.linalg.norm(axis) == 0:
            raise ValueError("zero-length C–X bond")
        cosine = np.divide(
            to_oxygen @ axis,
            length * np.linalg.norm(axis),
            out=np.ones_like(length),
            where=length > 0,
        )
        angle = np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0)))
        fraction = np.clip((angle - 135.0) / 30.0, 0.0, 1.0)
        angular = fraction * fraction * (3.0 - 2.0 * fraction)
        radial = np.exp(-0.5 * ((length - center) / WIDTH) ** 2)
        distance_only -= weight * radial.max()
        directional -= weight * (radial * angular).max()
    return distance_only, directional


def models(poses):
    """Recover the complete heavy-atom graph and correspondence from each model."""
    raw = poses.read_text().strip()
    if not raw.startswith("MODEL ") or not raw.endswith("ENDMDL"):
        raise ValueError(f"expected complete multi-model PDBQT: {poses}")
    for block in raw.split("ENDMDL")[:-1]:
        lines = block.strip().splitlines()
        if not lines[0].startswith("MODEL "):
            raise ValueError(f"missing MODEL record in {poses}")
        model = int(lines[0].split()[1])
        (score_line,) = [l for l in lines if l.startswith("REMARK VINA RESULT:")]
        (smiles_line,) = [
            l
            for l in lines
            if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX ")
        ]
        smiles = smiles_line[len("REMARK SMILES ") :]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"invalid Meeko graph in model {model}")
        types = {
            int(l[6:11]): l.split()[-1]
            for l in lines
            if l.startswith(("ATOM  ", "HETATM"))
        }
        coordinates = {
            int(l[6:11]): pdb_coordinates(l)
            for l in lines
            if l.startswith(("ATOM  ", "HETATM"))
        }
        mapping = {}
        for line in lines:
            if line.startswith("REMARK SMILES IDX "):
                fields = [int(n) for n in line.split()[3:]]
                if len(fields) % 2:
                    raise ValueError("odd REMARK SMILES IDX row")
                for smi, pdbqt in zip(fields[::2], fields[1::2]):
                    if smi - 1 in mapping or pdbqt in mapping.values():
                        raise ValueError("duplicate REMARK SMILES IDX correspondence")
                    mapping[smi - 1] = pdbqt
        if set(mapping) != set(range(mol.GetNumAtoms())):
            raise ValueError("incomplete REMARK SMILES IDX mapping")
        expected = {"A": "C", "NA": "N", "OA": "O", "SA": "S"}
        if any(
            expected.get(types[pdbqt], types[pdbqt])
            != mol.GetAtomWithIdx(i).GetSymbol()
            for i, pdbqt in mapping.items()
        ):
            raise ValueError("heavy atom types disagree with REMARK SMILES IDX mapping")
        yield (
            model,
            float(score_line.split()[3]),
            mol,
            np.array([coordinates[mapping[i]] for i in range(mol.GetNumAtoms())]),
        )


def rank_pool(receptor, poses):
    acceptors = carbonyls(receptor)
    rows = []
    for model, score, mol, xyz in models(poses):
        radial, directed = corrections(mol, xyz, acceptors)
        rows.append(
            dict(
                model=model,
                native=score,
                radial_delta=radial,
                directional_delta=directed,
                radial=score + radial,
                directional=score + directed,
            )
        )
    if not rows:
        raise ValueError("empty pose pool")
    for arm in ("native", "radial", "directional"):
        ordering = sorted(rows, key=lambda r: (r[arm], r["model"]))
        for rank, row in enumerate(ordering, 1):
            row[arm + "_rank"] = rank
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receptor", type=Path)
    parser.add_argument("poses", type=Path)
    args = parser.parse_args()
    rows = rank_pool(args.receptor, args.poses)
    w = csv.DictWriter(sys.stdout, fieldnames=rows[0])
    w.writeheader()
    w.writerows(rows)
