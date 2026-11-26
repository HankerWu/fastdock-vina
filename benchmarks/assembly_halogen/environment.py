"""Measure crystal packing contacts and a literature positive geometry in raw depositions."""

import argparse
import csv
import json
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem


def packing(folder):
    structure = gemmi.read_structure(str(folder / "deposition.pdb"))
    neighbors = gemmi.NeighborSearch(structure[0], structure.cell, 5.0).populate(
        include_h=False
    )
    ligand = Chem.SDMolSupplier(str(folder / "crystal.sdf"), removeHs=True)[0]
    observed = set()
    nearest = 100.0
    for index, xyz in enumerate(ligand.GetConformer().GetPositions()):
        point = gemmi.Position(*xyz)
        for mark in neighbors.find_atoms(point, radius=4):
            cra = mark.to_cra(structure[0])
            if cra.residue.het_flag != "A" or cra.atom.altloc not in ("\x00", "A"):
                continue
            here = structure.cell.find_nearest_pbc_position(
                point, cra.atom.pos, mark.image_idx
            )
            if here.dist(cra.atom.pos) < 1e-4:
                continue
            d = here.dist(point)
            if d <= 4:
                nearest = min(nearest, d)
                observed.add(index)
    return {
        "case": folder.name,
        "ligand_atoms_with_symmetry_protein_within_4A": len(observed),
        "minimum_image_distance_A": nearest if observed else "",
    }


def positive(path):
    pdb = gemmi.read_structure(str(path))[0]
    ligand = next(r for r in pdb["A"] if r.name == "1CX")
    x = next(a for a in ligand if a.element.name == "I")
    c = min(
        (a for a in ligand if a.element.name == "C"), key=lambda a: a.pos.dist(x.pos)
    )
    acceptor = min(
        ((a, r) for r in pdb["A"] if r.het_flag == "A" for a in r if a.name == "O"),
        key=lambda pair: pair[0].pos.dist(x.pos),
    )
    o, r = acceptor
    u = np.array(list(c.pos)) - np.array(list(x.pos))
    v = np.array(list(o.pos)) - np.array(list(x.pos))
    angle = np.degrees(np.arctan2(np.linalg.norm(np.cross(u, v)), np.dot(u, v)))
    return {
        "PDB": "3DY7",
        "ligand": ligand.name,
        "C": c.name,
        "I": x.name,
        "acceptor": f"{r.name} A {r.seqid.num} {o.name}",
        "distance_A": x.pos.dist(o.pos),
        "angle_deg": angle,
        "correction": -0.40
        * np.exp(-0.5 * ((x.pos.dist(o.pos) - 3.50) / 0.35) ** 2)
        * min(1, max(0, (angle - 135) / 30)) ** 2
        * (3 - 2 * min(1, max(0, (angle - 135) / 30))),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("positive", type=Path)
    p.add_argument("results", type=Path)
    a = p.parse_args()
    rows = [packing(f) for f in sorted(a.work.iterdir())]
    dest = a.results / "packing.csv"
    with dest.open("w") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)
    (a.results / "positive.json").write_text(
        json.dumps(positive(a.positive), indent=2) + "\n"
    )
