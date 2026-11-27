"""Fixed-coordinate, bidirectional aliphatic CH--aromatic carbon scoring."""

import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from rank_halogen_pool import AMINO_ACIDS, pdb_coordinates
from score_physical_pool import PROTEIN_RINGS, serial_index_pairs

# Eq. 4, Xiao & Woods, JCTC 2023, DOI:10.1021/acs.jctc.3c00300.
DEPTH, CENTER, WIDTH, WEIGHT = 0.26, 4.49, 0.75, 0.3
CUTOFF = 8.0
ATOMS = ("ATOM  ", "HETATM")

# Canonical protein residues are saturated except for the peptide/side-chain
# carbonyls, guanidinium, and rings in PROTEIN_RINGS.
SP2_SIDECHAINS = {
    "ARG": {"CZ"},
    "ASN": {"CG"},
    "ASP": {"CG"},
    "ASH": {"CG"},
    "GLN": {"CD"},
    "GLU": {"CD"},
    "GLH": {"CD"},
}
RESIDUES = AMINO_ACIDS | set(PROTEIN_RINGS)
RING_NAMES = {res: set().union(*rings) for res, rings in PROTEIN_RINGS.items()}


def kernel(d):
    """Raw, attractive pair contribution in kcal/mol at distance d in Å."""
    return np.where(
        d < CUTOFF, -DEPTH * np.exp(-0.5 * ((d - CENTER) / WIDTH) ** 2), 0.0
    )


def _graph_record(lines):
    return tuple(l for l in lines if l.startswith("REMARK SMILES "))


class CHPi:
    def __init__(self, receptor_text, ligand_template_text):
        r_arom, r_ch = [], []
        for line in receptor_text.splitlines():
            if not line.startswith(ATOMS):
                continue
            name, residue = line[12:16].strip(), line[17:20].strip()
            if residue not in RESIDUES or line.split()[-1] not in ("C", "A"):
                continue
            if name in RING_NAMES.get(residue, ()):
                r_arom.append(pdb_coordinates(line))
            elif name != "C" and name not in SP2_SIDECHAINS.get(residue, ()):
                r_ch.append(pdb_coordinates(line))
        self.r_arom = np.asarray(r_arom, dtype=float).reshape(-1, 3)
        self.r_ch = np.asarray(r_ch, dtype=float).reshape(-1, 3)
        if not (np.isfinite(self.r_arom).all() and np.isfinite(self.r_ch).all()):
            raise ValueError("nonfinite protein coordinates")

        lines = ligand_template_text.splitlines()
        self.graph = _graph_record(lines)
        (smiles,) = [
            l.split()[2] for l in self.graph if not l.startswith("REMARK SMILES IDX ")
        ]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("invalid ligand SMILES")
        pairs = serial_index_pairs(lines, "REMARK SMILES IDX ")
        serial = dict(pairs)
        n = mol.GetNumAtoms()
        if (
            len(pairs) != n
            or set(serial) != set(range(1, n + 1))
            or len(set(serial.values())) != n
        ):
            raise ValueError("ligand requires a complete one-to-one chemical graph map")
        self.l_arom, self.l_ch = [], []
        types = {int(l[6:11]): l.split()[-1] for l in lines if l.startswith(ATOMS)}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            i = serial[atom.GetIdx() + 1]
            if types.get(i) not in ("C", "A"):
                raise ValueError(
                    f"mapped carbon {i} has unsupported atom type {types.get(i)}"
                )
            if atom.GetIsAromatic():
                self.l_arom.append(i)
            elif (
                atom.GetHybridization() == Chem.HybridizationType.SP3
                and atom.GetTotalNumHs() > 0
            ):
                self.l_ch.append(i)

    def components(self, pose_text):
        """Ligand CH→protein aromatic, then protein CH→ligand aromatic; unweighted."""
        lines = pose_text.splitlines()
        graph = _graph_record(lines)
        if graph and graph != self.graph:
            raise ValueError(
                "pose topology or serial correspondence differs from template"
            )
        records = [line for line in lines if line.startswith(ATOMS)]
        xyz = {int(line[6:11]): pdb_coordinates(line) for line in records}
        if len(xyz) != len(records):
            raise ValueError("duplicate ligand serial or multiple models")

        def energy(ids, partner):
            points = np.asarray([xyz[i] for i in ids], dtype=float).reshape(-1, 3)
            if not np.isfinite(points).all():
                raise ValueError("nonfinite ligand coordinates")
            return float(kernel(cdist(points, partner)).sum())

        return energy(self.l_ch, self.r_arom), energy(self.l_arom, self.r_ch)
