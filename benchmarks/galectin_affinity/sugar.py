"""Chemical correspondence of the cyclic donor and query sugar heavy atoms."""

import numpy as np
from rdkit import Chem

DONORS = {
    "1GZW": ("D", "GAL"),
    "4Q26": ("C", "GAL"),
    "1KJL": ("B", "GAL"),
    "3T1L": ("A", "MQT"),
}
ATOM_NAMES = ("C1", "C2", "C3", "C4", "C5", "O5")
CORE = Chem.MolFromSmarts("[C;R]1([S])[C;R][C;R][C;R][C;R][O;R]1")


def xyz(line):
    return [float(line[i : i + 8]) for i in (30, 38, 46)]


def crystal_core(work, pdb):
    chain, comp = DONORS[pdb]
    rows = (work / "raw" / f"{pdb}.pdb").read_text().splitlines()
    atoms = {
        r[12:16].strip(): xyz(r)
        for r in rows
        if r.startswith("HETATM")
        and r[21] == chain
        and r[17:20].strip() == comp
        and r[16] in " A"
    }
    return np.array([atoms[a] for a in ATOM_NAMES])


def ligand_core(work, compound):
    rows = (work / "prepared" / f"{compound}.pdbqt").read_text().splitlines()
    (smiles,) = [
        s.split()[2]
        for s in rows
        if s.startswith("REMARK SMILES ") and not s.startswith("REMARK SMILES IDX")
    ]
    mol = Chem.MolFromSmiles(smiles)
    (match,) = mol.GetSubstructMatches(CORE)
    serials = {}
    for r in rows:
        if r.startswith("REMARK SMILES IDX"):
            values = list(map(int, r.split()[3:]))
            serials.update(zip(values[::2], values[1::2]))
    return [serials[i + 1] for i in (match[0], *match[2:])]
