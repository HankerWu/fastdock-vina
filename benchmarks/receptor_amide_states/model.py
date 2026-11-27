"""Fixed-coordinate, two-state amide identity replacement for the Vina score."""

import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from benchmarks.frozen_core.hydroxyl_partition.scorer import (
    ACCEPTORS,
    RT,
    atoms,
    ligand_sites,
    receptor_hydrogens,
)

RADII = {
    "A": 1.9,
    "C": 1.9,
    "N": 1.8,
    "NA": 1.8,
    "O": 1.7,
    "OA": 1.7,
    "S": 2.0,
    "SA": 2.0,
    "P": 2.1,
    "F": 1.5,
    "Cl": 1.8,
    "Br": 2.0,
    "I": 2.2,
    "Si": 2.2,
    "At": 2.3,
}
METALS = {"Mg", "Mn", "Zn", "Ca", "Fe", "Cu", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"}
RADII.update(dict.fromkeys(METALS, 1.2))
RADII.update(dict.fromkeys(("CG0", "CG1", "CG2", "CG3"), 1.9))
IGNORED = {"HD", "H", "G0", "G1", "G2", "G3", "W"}
# Within the same residue remove a terminus's 1-2 to 1-4 neighbors. Backbone
# amides can H-bond back to their own Asn/Gln side chain and are retained.
TERMINI = {
    "ASN": ("OD1", "ND2", {"CG", "CB", "CA"}),
    "GLN": ("OE1", "NE2", {"CD", "CG", "CB"}),
}


def sigma(kind):
    return RADII[{"Se": "S"}.get(kind, kind)]


def pair_energy(distances, terminal_radius, neighbor_radii, can_hbond):
    """Nonhydrophobic Vina pair potential, before native ligand normalization."""
    d = distances - terminal_radius - neighbor_radii
    e = (
        -0.035579 * np.exp(-((d / 0.5) ** 2))
        - 0.005156 * np.exp(-(((d - 3.0) / 2.0) ** 2))
        + 0.840245 * np.minimum(d, 0.0) ** 2
        - 0.587439 * can_hbond * np.clip(-d / 0.7, 0.0, 1.0)
    )
    return np.where(distances < 8.0, e, 0.0)


def swap_delta(distances, radii, donor, acceptor):
    """Difference flipped minus deposited, per external atom. O position first."""
    at_o, at_n = distances
    return (
        pair_energy(at_o, 1.8, radii, acceptor) - pair_energy(at_o, 1.7, radii, donor)
    ) + (pair_energy(at_n, 1.7, radii, donor) - pair_energy(at_n, 1.8, radii, acceptor))


def marginal(apo_delta, ligand_delta):
    return RT * (
        np.logaddexp(0.0, -apo_delta / RT)
        - np.logaddexp(0.0, -(apo_delta + ligand_delta) / RT)
    )


@dataclass
class Site:
    residue: str
    endpoints: np.ndarray
    apo_delta: float


class AmideReceptor:
    def __init__(self, pdbqt):
        parsed = atoms(pdbqt)
        residues = defaultdict(list)
        for serial, atom in parsed.items():
            residues[atom.residue].append(serial)
        owner = receptor_hydrogens(parsed, residues)
        ids = [serial for serial, atom in parsed.items() if atom.kind not in IGNORED]
        heavy = [parsed[serial] for serial in ids]
        xyz = np.array([a.xyz for a in heavy])
        tree = cKDTree(xyz)
        radii = np.array([sigma(a.kind) for a in heavy])
        donor = np.array(
            [
                a.kind in METALS
                or (a.kind in {"NA", "N", "OA", "O"} and bool(owner.get(i)))
                for i, a in zip(ids, heavy)
            ]
        )
        acceptor = np.array([a.kind in ACCEPTORS for a in heavy])
        self.sites = []
        self.excluded = []
        for residue, group in residues.items():
            if residue[:3] not in TERMINI:
                continue
            o_name, n_name, excluded_names = TERMINI[residue[:3]]
            names = {parsed[i].name: i for i in group}
            if (
                o_name not in names
                or n_name not in names
                or parsed[names[o_name]].kind != "OA"
                or parsed[names[n_name]].kind != "N"
                or len(owner.get(names[n_name], ())) != 2
            ):
                self.excluded.append(residue)
                continue
            endpoints = np.array([parsed[names[name]].xyz for name in (o_name, n_name)])
            nearby = sorted(set().union(*tree.query_ball_point(endpoints, 8.0)))
            selected = [
                i
                for i in nearby
                if heavy[i].residue != residue
                or heavy[i].name not in excluded_names | {o_name, n_name}
            ]
            ds = np.linalg.norm(endpoints[:, None, :] - xyz[selected], axis=2)
            delta = swap_delta(
                ds, radii[selected], donor[selected], acceptor[selected]
            ).sum()
            self.sites.append(Site(residue, endpoints, float(delta)))
        self.endpoints = np.array([s.endpoints for s in self.sites]).reshape(-1, 2, 3)
        self.state_energy = np.array([s.apo_delta for s in self.sites])

    def score(self, ligand_text, block, *, detail=False):
        template, owners = ligand_sites(ligand_text)
        pose = atoms(block)
        if pose.keys() != template.keys() or any(
            a.kind != template[i].kind for i, a in pose.items()
        ):
            raise ValueError("Pose does not match prepared ligand serials/types")
        ids = [i for i, a in pose.items() if a.kind not in IGNORED]
        kinds = [pose[i].kind for i in ids]
        xyz = np.array([pose[i].xyz for i in ids]).reshape(-1, 3)
        distances = np.linalg.norm(self.endpoints[:, :, None, :] - xyz, axis=3)
        active = np.flatnonzero((distances < 8.0).any(axis=(1, 2)))
        radii = np.array([sigma(k) for k in kinds])
        donors = np.array(
            [
                k in METALS or (k in {"NA", "N", "OA", "O"} and bool(owners.get(i)))
                for i, k in zip(ids, kinds)
            ],
            dtype=bool,
        )
        acceptors = np.array([k in ACCEPTORS for k in kinds], dtype=bool)
        delta = swap_delta(
            distances[active].transpose(1, 0, 2), radii, donors, acceptors
        ).sum(axis=1)
        changes = marginal(self.state_energy[active], delta)
        terms = {
            "ensemble": float(changes.sum()),
            "minimum": float(np.minimum(delta, 0.0).sum()),
            "uniform": float(marginal(0.0, delta).sum()),
            "flip_all": float(delta.sum()),
            "active_sites": len(active),
        }
        if detail:
            terms["sites"] = [
                (self.sites[i].residue, self.sites[i].apo_delta, float(dl), float(c))
                for i, dl, c in zip(active, delta, changes)
            ]
        return terms
