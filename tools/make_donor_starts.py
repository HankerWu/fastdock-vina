"""Place an independently prepared ligand using rotatable dihedrals from a bound analogue."""

import argparse
import itertools
from bisect import insort
from pathlib import Path

import numpy as np
from make_scaffold_anchors import all_correspondences
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import rdMolTransforms


def make_starts(query, donor, copy_torsions=True, limit=4):
    """Return (mapped RMSD, atom assignment, positioned all-H query) in donor frame."""
    if limit < 1 or not all(
        m and m.GetNumConformers() == 1 and m.GetConformer().Is3D()
        for m in (query, donor)
    ):
        raise ValueError(
            "Expected two single-conformer 3D molecules and a positive limit"
        )
    core = Chem.RemoveHs(query)
    query = Chem.AddHs(core, addCoords=True)
    donor = Chem.RemoveHs(donor)
    points = donor.GetConformer().GetPositions()
    mappings = all_correspondences(core, donor)
    (setup,) = MoleculePreparation().prepare(query)
    tree, ok, error = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(error)
    ids = {}
    for l in tree.splitlines():
        if l.startswith("REMARK INDEX MAP"):
            v = list(map(int, l.split()[3:]))
            ids.update((s, i - 1) for i, s in zip(v[::2], v[1::2]))
    axes = [
        tuple(ids[int(i)] for i in l.split()[1:])
        for l in tree.splitlines()
        if l.startswith("BRANCH ")
    ]
    candidates = []
    for mapping in mappings:
        map_dict = dict(mapping)
        mol = Chem.Mol(query)
        conf = mol.GetConformer()
        if copy_torsions:
            for a, b in axes:
                if a not in map_dict or b not in map_dict:
                    continue
                left = sorted(
                    n.GetIdx()
                    for n in core.GetAtomWithIdx(a).GetNeighbors()
                    if n.GetIdx() != b and n.GetIdx() in map_dict
                )
                right = sorted(
                    n.GetIdx()
                    for n in core.GetAtomWithIdx(b).GetNeighbors()
                    if n.GetIdx() != a and n.GetIdx() in map_dict
                )
                if not left or not right:
                    continue
                q = (left[0], a, b, right[0])
                d = [map_dict[i] for i in q]
                if any(
                    donor.GetBondBetweenAtoms(i, j) is None
                    for i, j in itertools.pairwise(d)
                ):
                    continue
                angle = rdMolTransforms.GetDihedralDeg(donor.GetConformer(), *d)
                if not np.isfinite(angle):
                    raise ValueError("Degenerate donor dihedral")
                rdMolTransforms.SetDihedralDeg(conf, *q, angle)
        qi, di = zip(*mapping)
        xyz = conf.GetPositions()
        x, y = xyz[list(qi)], points[list(di)]
        u, _, vt = np.linalg.svd((x - x.mean(0)).T @ (y - y.mean(0)))
        u[:, -1] *= np.sign(np.linalg.det(u @ vt))
        placed = (xyz - x.mean(0)) @ u @ vt + y.mean(0)
        error = float(np.sqrt(np.mean(np.sum((placed[list(qi)] - y) ** 2, axis=1))))
        for i, pt in enumerate(placed):
            conf.SetAtomPosition(i, pt)
        insort(candidates, (error, mapping, mol), key=lambda c: (c[0], c[1]))
        if len(candidates) > limit:
            candidates.pop()
    return candidates


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "query", type=Path, help="SDF with an independent starting conformation"
    )
    p.add_argument(
        "donor", type=Path, help="SDF of a different bound ligand in the receptor frame"
    )
    p.add_argument("output", type=Path)
    p.add_argument(
        "--rigid-only", action="store_true", help="Alignment-only experimental control"
    )
    p.add_argument("--limit", type=int, default=4)
    a = p.parse_args()
    records = [
        list(Chem.SDMolSupplier(str(f), removeHs=False)) for f in (a.query, a.donor)
    ]
    if any(len(r) != 1 or r[0] is None for r in records):
        p.error("Each file must contain exactly one valid molecule")
    if a.output.resolve() in (a.query.resolve(), a.donor.resolve()):
        p.error("Output must differ from both inputs")
    starts = make_starts(records[0][0], records[1][0], not a.rigid_only, a.limit)
    with Chem.SDWriter(str(a.output)) as w:
        for error, mapping, mol in starts:
            mol.SetDoubleProp("donor_core_rmsd", error)
            mol.SetIntProp("matched_atoms", len(mapping))
            w.write(mol)


if __name__ == "__main__":
    main()
