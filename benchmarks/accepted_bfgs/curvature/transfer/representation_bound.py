"""Necessary RMSD bounds for the prepared ligand, with or without amide rotations."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.optimize import linprog

from evaluate import PDBS, write_csv

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "receptor_linkage/aurora"))
from pose_probe import get_tree, register


def rigid_blocks(heavy, rotors):
    """Maximal heavy-atom subsets carried together by every rotation.

    The two axis endpoints belong to both sides since rotating them makes
    no change. Intersecting the possible sides for every independent rotor
    yields sets that are invariant under a sequence of such rotations.
    """
    atoms = set(heavy)
    blocks = [atoms]
    for a, b, subtree in rotors:
        inside = set(subtree)
        sides = (inside | {a, b}, (atoms - inside) | {a, b})
        candidates = {frozenset(block & side) for block in blocks for side in sides}
        blocks = [s for s in candidates if s and not any(s < t for t in candidates)]
    return sorted(sorted(heavy.index(i) for i in block) for block in blocks)


def bound(points, target, blocks):
    # A block can receive its own optimal translation/rotation in this relaxation.
    matrix = np.zeros((len(blocks), len(points)))
    residual = []
    for k, ix in enumerate(blocks):
        matrix[k, ix] = 1
        aligned = register(points, ix, target[ix])
        residual.append(np.square(aligned[ix] - target[ix]).sum())
    # The squared displacement of each atom must contribute at least the block
    # residual in every block containing it. The LP handles shared axis atoms.
    sol = linprog(np.ones(len(points)), A_ub=-matrix, b_ub=-np.array(residual),
                  bounds=(0, None), method="highs")
    if not sol.success:
        raise RuntimeError(sol.message)
    weights = np.maximum(0, -sol.ineqlin.marginals)
    weights /= max(1, np.max(matrix.T @ weights))
    dual = weights @ residual
    if abs(dual - sol.fun) > 1e-7:
        raise ArithmeticError("Primal and independently feasible dual bounds differ")
    return np.sqrt(dual / len(points))


def evaluate(work, output):
    rows = []
    for pdb in PDBS:
        meta = json.loads((work / pdb / "metadata.json").read_text())
        _, serials, sdf_ids, coords, arms = get_tree(work / pdb / "ligand.pdbqt")
        heavy = [serials.index(s) for s in meta["serial"]]
        mol = Chem.MolFromSmiles(meta["smiles"])
        amides = {frozenset((c,n)) for c,_,n in mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)N"))}
        free = [arm for arm in arms if frozenset(sdf_ids[i] for i in arm[:2]) not in amides]
        for label, allowed in (("declared", arms), ("all_amides_locked", free)):
            blocks = rigid_blocks(heavy, allowed)
            values = [bound(coords[heavy], np.array(target), blocks) for target in meta["crystal"]]
            rows.append(dict(query=pdb, domain=label, heavy_atoms=len(heavy), axes=len(allowed),
                             blocks=len(blocks), symmetry_mappings=len(values), lower_bound_A=min(values)))
    write_csv(output, rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    evaluate(a.work, a.output)
