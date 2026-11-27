"""Empirical receptor O-H orientational marginal on fixed PDBQT coordinates."""
from collections import defaultdict
from functools import lru_cache
from typing import NamedTuple

import numpy as np
from scipy.spatial import cKDTree
from scipy.special import logsumexp

RT = 0.00198720425864083 * 298.15  # kcal/mol
ACCEPTORS = {"NA", "OA"}
RADII = {"NA": 1.8, "N": 1.8, "OA": 1.7, "O": 1.7}
AXES = {("SER", "OG"): "CB", ("THR", "OG1"): "CB", ("TYR", "OH"): "CZ"}


class Atom(NamedTuple):
    kind: str
    residue: str
    name: str
    xyz: np.ndarray


def atoms(text):
    parsed = {}
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            if serial in parsed:
                raise ValueError("Duplicate PDBQT serial; provide exactly one pose")
            parsed[serial] = Atom(line[77:].strip(), line[17:27], line[12:16].strip(),
                                  np.array([float(line[i:i+8]) for i in (30, 38, 46)]))
    return parsed


def magnitude(left, right):
    """Magnitude of Vina's native favorable N/O hbond term in raw score units."""
    overlap = RADII[left.kind] + RADII[right.kind] - np.linalg.norm(left.xyz - right.xyz)
    return .587439 * np.clip(overlap / .7, 0, 1)


def quality(donor, hydrogen, acceptor):
    dh = hydrogen - donor
    ha = acceptor - hydrogen
    norm = np.linalg.norm(dh, axis=-1) * np.linalg.norm(ha, axis=-1)
    dot = (dh * ha).sum(axis=-1)
    return np.clip(np.divide(dot, norm, out=np.zeros_like(dot), where=norm > 1e-12), 0, 1) ** 2


def circular_positions(carbon, oxygen, hydrogen, n, phase):
    axis = oxygen - carbon
    axis /= np.linalg.norm(axis)
    oh = hydrogen - oxygen
    parallel = axis * np.dot(axis, oh)
    u = oh - parallel
    v = np.cross(axis, u)
    phi = (np.arange(n) + phase) * (2 * np.pi / n)
    return oxygen + parallel + np.cos(phi)[:, None] * u + np.sin(phi)[:, None] * v


def marginal(penalties, prior_energies):
    """-RT log <exp(-penalty/RT)> under the specified normalized finite prior."""
    log_prior = -prior_energies / RT
    return RT * (logsumexp(log_prior) - logsumexp(log_prior - penalties / RT))


def receptor_hydrogens(parsed, residues):
    owners = defaultdict(list)
    for i, atom in parsed.items():
        if atom.kind == "HD":
            heavy = (j for j in residues[atom.residue] if parsed[j].kind not in {"H", "HD"})
            near = sorted((np.linalg.norm(atom.xyz - parsed[j].xyz), j) for j in heavy)
            if not near or near[0][0] > 1.5 or (len(near)>1 and near[1][0] - near[0][0]<.1):
                raise ValueError(f"Ambiguous or missing covalent receptor parent for HD {i}")
            owners[near[0][1]].append(i)
    return owners


class Receptor:
    def __init__(self, receptor_text, n=72, phase=0.):
        if n < 2:
            raise ValueError("At least two quadrature nodes required")
        self.atoms = atoms(receptor_text)
        residues = defaultdict(list)
        for i, atom in self.atoms.items():
            residues[atom.residue].append(i)
        self.hydrogens = receptor_hydrogens(self.atoms, residues)
        self.ids = [i for i, a in self.atoms.items() if a.kind in RADII]
        self.tree = cKDTree(np.array([self.atoms[i].xyz for i in self.ids]).reshape(-1, 3))
        self.rotors = {}
        for i, hs in self.hydrogens.items():
            oxygen = self.atoms[i]
            base = AXES.get((oxygen.residue[:3], oxygen.name))
            if base is None:
                continue
            carbons = [self.atoms[j] for j in residues[oxygen.residue] if self.atoms[j].name == base]
            if (len(hs) != 1 or len(carbons) != 1 or carbons[0].kind not in ("C", "A")
                    or oxygen.kind not in ("O", "OA") or not 1.<np.linalg.norm(oxygen.xyz-carbons[0].xyz)<1.8):
                raise ValueError(f"Incompatible hydroxyl carbon or hydrogen topology at {oxygen.residue}")
            hx = circular_positions(carbons[0].xyz, oxygen.xyz, self.atoms[hs[0]].xyz, n, phase)
            prior = np.zeros(n)
            for neighbor in self.tree.query_ball_point(oxygen.xyz, 3.6):
                acceptor = self.atoms[self.ids[neighbor]]
                if acceptor.residue != oxygen.residue and acceptor.kind in ACCEPTORS:
                    prior -= magnitude(oxygen, acceptor) * quality(oxygen.xyz, hx, acceptor.xyz)
            self.rotors[i] = (hx, prior)


@lru_cache(maxsize=128)
def ligand_sites(ligand_text):
    topology = atoms(ligand_text)
    mapped, pairs = {}, []
    for line in ligand_text.splitlines():
        for prefix, target in (("REMARK SMILES IDX ", mapped), ("REMARK H PARENT ", pairs)):
            if line.startswith(prefix):
                values = list(map(int, line[len(prefix):].split()))
                data = list(zip(values[::2], values[1::2], strict=True))
                if target is mapped:
                    mapped.update(data)
                else:
                    pairs.extend(data)
    owners = defaultdict(list)
    for parent, serial in pairs:
        owners[mapped[parent]].append(serial)
    if sorted(h for hs in owners.values() for h in hs) != sorted(i for i,a in topology.items() if a.kind == "HD"):
        raise ValueError("Ligand needs complete Meeko REMARK H PARENT / SMILES IDX records")
    return topology, owners


def score_pose(receptor, ligand_text, pose_block):
    topology, lh = ligand_sites(ligand_text)
    xyz = atoms(pose_block)
    if xyz.keys() != topology.keys() or any(a.kind != topology[i].kind for i, a in xyz.items()):
        raise ValueError("Pose serials or atom types differ from the prepared query")
    nonrotor = radial = fixed_rotor = 0.
    contacts = mobile_contacts = 0
    fields = {}
    for i, ligand in xyz.items():
        if ligand.kind not in RADII:
            continue
        for offset in receptor.tree.query_ball_point(ligand.xyz, 3.6):
            j = receptor.ids[offset]
            protein = receptor.atoms[j]
            weight = magnitude(ligand, protein)
            lig_h = lh[i] if protein.kind in ACCEPTORS else ()
            rec_h = receptor.hydrogens[j] if ligand.kind in ACCEPTORS else ()
            if weight <= 0 or not (lig_h or rec_h):
                continue
            contacts += 1
            radial -= weight
            forward = max((quality(ligand.xyz, xyz[h].xyz, protein.xyz) for h in lig_h), default=0.)
            back = max((quality(protein.xyz, receptor.atoms[h].xyz, ligand.xyz) for h in rec_h), default=0.)
            fixed = weight * (1 - max(forward, back))
            if rec_h and j in receptor.rotors:
                mobile_contacts += 1
                fixed_rotor += fixed
                positions, _ = receptor.rotors[j]
                penalty = weight * (1 - np.maximum(forward, quality(protein.xyz, positions, ligand.xyz)))
                if j in fields:
                    fields[j] += penalty
                else:
                    fields[j] = penalty
            else:
                nonrotor += fixed
    minimum = uniform = ensemble = 0.
    for j, penalty in fields.items():
        _, prior = receptor.rotors[j]
        minimum += penalty.min()
        uniform += marginal(penalty, np.zeros_like(prior))
        ensemble += marginal(penalty, prior)
    return dict(contacts=contacts, rotors=len(fields), mobile_contacts=mobile_contacts, raw_radial=radial,
                raw_fixed=nonrotor+fixed_rotor, raw_minimum=nonrotor+minimum,
                raw_uniform=nonrotor+uniform, raw_ensemble=nonrotor+ensemble, raw_rotor_only=ensemble,
                raw_nonrotor=nonrotor)


if __name__ == "__main__":
    import argparse
    import csv
    from pathlib import Path
    import sys

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receptor", type=Path)
    parser.add_argument("ligand", type=Path)
    parser.add_argument("pool", type=Path, help="Scored multi-model Vina PDBQT")
    parser.add_argument("--scale", type=float, required=True,
                        help="Use the native scoring binary's conf_independent(ligand,1) value")
    args = parser.parse_args()
    if not 0 < args.scale <= 1:
        parser.error("Expected a positive native Vina torsional reporting scale at most one")
    rec = Receptor(args.receptor.read_text())
    lig = args.ligand.read_text()
    rows = []
    for rank, block in enumerate((b for b in args.pool.read_text().split("ENDMDL") if b.strip()),1):
        score = float(next(l for l in block.splitlines() if l.startswith("REMARK VINA RESULT:")).split()[3])
        terms = score_pose(rec, lig, block)
        rows.append(dict(rank=rank, native=score, predicted=score+args.scale*terms["raw_ensemble"], **terms))
    if not rows:
        parser.error("The pose pool is empty")
    writer = csv.DictWriter(sys.stdout,fieldnames=rows[0])
    writer.writeheader()
    writer.writerows(rows)
