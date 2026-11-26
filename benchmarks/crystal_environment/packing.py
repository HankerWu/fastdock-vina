"""Protein crystal neighbours and a Vina repulsive-energy diagnostic."""

import itertools
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree

RADII = dict(C=1.9, N=1.8, O=1.7, S=2.0, P=2.1, F=1.5, CL=1.8, BR=2.0, I=2.2, SE=2.0)
WEIGHT = 0.840245


@dataclass(frozen=True)
class Atom:
    name: str
    residue: str
    chain: str
    number: str
    alt: str
    element: str
    occ: float
    b: float
    xyz: tuple
    polymer: bool


def atoms(path):
    values = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        e = line[76:78].strip().upper()
        if e in ("H", "D") or float(line[54:60]) <= 0:
            continue
        values.append(
            Atom(
                line[12:16].strip(),
                line[17:20].strip(),
                line[21],
                line[22:27].strip(),
                line[16].strip(),
                e,
                float(line[54:60]),
                float(line[60:66]),
                tuple(float(line[i : i + 8]) for i in (30, 38, 46)),
                line.startswith("ATOM  ") or line[17:20] == "MSE",
            )
        )
    return values


def crystal_neighbors(path, host, query, padding=18.0):
    """Enumerate all ASU copies possibly within padding of the query Cartesian box."""
    st = gemmi.read_structure(str(path))
    cell = np.asarray(st.cell.orth.mat, dtype=float)
    frac = np.linalg.inv(cell)
    seq = [a for a in atoms(path) if a.polymer]
    xyz = np.array([a.xyz for a in seq])
    xyz_f = xyz @ frac.T
    q_f = query @ frac.T
    # A Cartesian sphere of radius padding has projection padding*|b_i| in fractional axis i.
    margin = padding * np.linalg.norm(frac, axis=1)
    center_min, center_max = (
        np.min(query, axis=0) - padding,
        np.max(query, axis=0) + padding,
    )
    positions, records, images = [], [], []
    for op in st.find_spacegroup().operations():
        rot = np.array(op.rot) / op.DEN
        shift = np.array(op.tran) / op.DEN
        mapped = xyz_f @ rot.T + shift
        lo = np.ceil(q_f.min(0) - margin - mapped.max(0)).astype(int)
        hi = np.floor(q_f.max(0) + margin - mapped.min(0)).astype(int)
        for lattice in itertools.product(*(range(l, h + 1) for l, h in zip(lo, hi))):
            ident = op.triplet() == "x,y,z" and lattice == (0, 0, 0)
            cart = (mapped + lattice) @ cell.T
            near = ((cart >= center_min) & (cart <= center_max)).all(axis=1)
            for i in np.flatnonzero(near):
                if ident and seq[i].chain == host:
                    continue
                positions.append(cart[i])
                records.append(seq[i])
                images.append(f"{op.triplet()}|{','.join(map(str, lattice))}")
    return Environment(positions, records, images)


class Environment:
    def __init__(self, coords, records, images):
        self.coords = np.asarray(coords, dtype=float).reshape(-1, 3)
        self.records = records
        self.images = images
        self.tree = cKDTree(self.coords)
        self.radii = np.array([RADII[a.element] for a in records])
        self.occ = np.array([a.occ for a in records])

    def contacts(self, xyz, elements):
        distances, nearby = self.tree.query(xyz)
        i = int(np.argmin(distances))
        hit = self.records[nearby[i]] if len(self.coords) else None
        lig_tree = cKDTree(xyz)
        neighbors = lig_tree.sparse_distance_matrix(
            self.tree, 4.5, output_type="ndarray"
        )
        row, col, dist = (neighbors[k] for k in ("i", "j", "v"))
        sigma = np.array([RADII[e] for e in elements])[row] + self.radii[col]
        rep = WEIGHT * np.sum(self.occ[col] * np.square(np.maximum(sigma - dist, 0)))
        return dict(
            closest=float(distances[i]),
            n4=int((dist <= 4).sum()),
            nclash=int((dist < 2.3).sum()),
            repulsion=float(rep),
            ligand_nearest=i,
            nearest_host=f"{hit.chain}:{hit.residue}{hit.number}:{hit.name}"
            if hit
            else "",
            image=self.images[nearby[i]] if hit else "",
        )
