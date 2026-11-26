"""Create two-point or two-distance search anchors for an aryl amidinium and a protein carboxylate."""

import argparse
from pathlib import Path

import numpy as np
from rdkit import Chem
from make_scaffold_anchors import format_assignments
from score_physical_pool import serial_index_pairs

AMIDINE = Chem.MolFromSmarts("[c]-[C;X3](=[N+;H2])-[N;H2]")
ACID = Chem.MolFromSmarts("[C;X3](=O)-[O-]")
GAP = 2.90


def xyz(line):
    return [float(line[k : k + 8]) for k in (30, 38, 46)]


def site_points(c, oxygens):
    """Continue the carboxylate bisector toward a planar, bifurcated ligand site."""
    bisector = np.mean(oxygens, axis=0) - c
    length = np.linalg.norm(bisector)
    if length < 0.1:
        raise ValueError("Undefined carboxylate bisector")
    return oxygens + GAP * bisector / length


def carboxylate(path, chain, residue, insertion=""):
    """Use a prepared, single receptor geometry with unique carboxylate atoms."""
    records = {}
    residue_name = None
    for line in path.read_text().splitlines():
        if (
            not line.startswith(("ATOM  ", "HETATM"))
            or line[21] != chain
            or int(line[22:26]) != residue
            or line[26].strip() != insertion
        ):
            continue
        name = line[17:20]
        if name not in ("ASP", "GLU"):
            continue
        if residue_name and name != residue_name:
            raise ValueError("Ambiguous residue identity")
        residue_name = name
        atom = line[12:16].strip()
        if atom not in ("CG", "OD1", "OD2", "CD", "OE1", "OE2"):
            continue
        if line[16] not in (" ", "A"):
            continue
        if atom in records:
            raise ValueError(f"Ambiguous carboxylate atom {atom}")
        records[atom] = xyz(line)
    names = ("CG", "OD1", "OD2") if residue_name == "ASP" else ("CD", "OE1", "OE2")
    if any(name not in records for name in names):
        raise ValueError(
            "Missing carboxylate atoms for the specified residue and chain"
        )
    acid = np.asarray([records[name] for name in names])
    if (
        not np.isfinite(acid).all()
        or min(np.linalg.norm(acid[1:] - acid[0], axis=1)) < 0.5
    ):
        raise ValueError("Invalid receptor carboxylate coordinates")
    return acid


def ligand_sites(ligand):
    lines = ligand.read_text().splitlines()
    smiles = [
        line.split()[2]
        for line in lines
        if line.startswith("REMARK SMILES ")
        and not line.startswith("REMARK SMILES IDX")
    ]
    if len(smiles) != 1:
        raise ValueError(f"Expected a single Meeko SMILES record in {ligand}")
    mol = Chem.MolFromSmiles(smiles[0])
    if mol is None:
        raise ValueError(f"Unreadable Meeko SMILES in {ligand}")
    serials = {i - 1: s for i, s in serial_index_pairs(lines, "REMARK SMILES IDX ")}
    if len(serials) != mol.GetNumHeavyAtoms() or len(set(serials.values())) != len(
        serials
    ):
        raise ValueError(f"Incomplete atom identifiers in {ligand}")
    positions = {
        int(line[6:11]): np.array(xyz(line))
        for line in lines
        if line.startswith(("ATOM", "HETATM"))
    }
    matches = mol.GetSubstructMatches(AMIDINE)
    if not matches:
        raise ValueError(f"No protonated aryl primary amidine in {ligand}")
    for _, carbon, *nitrogens in matches:
        for n in nitrogens:
            length = np.linalg.norm(positions[serials[n]] - positions[serials[carbon]])
            if not 1.0 < length < 1.67:
                raise ValueError(f"Inconsistent amidine serials in {ligand}: {length}")
    return [tuple(serials[x] for x in match[1:]) for match in matches], [
        tuple(serials[x] for x in match) for match in mol.GetSubstructMatches(ACID)
    ]


def anchors(ligand, acid, geometry):
    """Return all admissible group assignments for one named anionic receptor site."""
    groups, _ = ligand_sites(ligand)
    pair, point = [], []
    for _, n1, n2 in groups:
        for a, b in ((0, 1), (1, 0)):
            pair.append(((n1, a), (n2, b)))
            point.extend((((n1, a),), ((n2, b),)))
    serials = {s: s for group in groups for s in group}
    if geometry in ("point", "pair"):
        return format_assignments(
            point if geometry == "point" else pair,
            serials,
            site_points(acid[0], acid[1:]),
        )
    if geometry != "shell":
        raise ValueError(f"Unknown geometry {geometry}")
    blocks = []
    for assignment in pair:
        blocks.append(
            "\n".join(
                f"{n} "
                + " ".join(f"{v:.8f}" for v in acid[1 + oi])
                + f" shell {GAP:.2f}"
                for n, oi in assignment
            )
        )
    return "\nOR\n".join(blocks) + "\n"


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "receptor",
        type=Path,
        help="Prepared PDB or PDBQT in the docking reference frame",
    )
    p.add_argument(
        "ligand",
        type=Path,
        help="Meeko PDBQT with declared protonated primary aryl amidinium",
    )
    p.add_argument("out", type=Path)
    p.add_argument("--chain", required=True)
    p.add_argument("--residue", required=True, type=int)
    p.add_argument("--insertion", default="")
    p.add_argument("--geometry", required=True, choices=("point", "pair", "shell"))
    a = p.parse_args()
    try:
        acid = carboxylate(a.receptor, a.chain, a.residue, a.insertion)
        a.out.write_text(anchors(a.ligand, acid, a.geometry))
    except (OSError, ValueError) as err:
        p.exit(2, f"{p.prog}: {err}\n")
