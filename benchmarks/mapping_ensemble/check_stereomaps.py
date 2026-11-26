"""Audit tetrahedral and alkene map orientation using RDKit and frozen real graphs."""

import argparse
import json
import random
import sys
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from make_scaffold_anchors import (
    all_correspondences,
    stereo_compatible,
    stereo_constraints,
)


def check_invariance(query, donor, rng):
    direct = set(all_correspondences(query, donor))
    swapped = {
        tuple(sorted((d, q) for q, d in pairs))
        for pairs in all_correspondences(donor, query)
    }
    assert swapped == direct, "Donor/query exchange changed valid maps"
    orders = [
        rng.sample(range(m.GetNumAtoms()), m.GetNumAtoms()) for m in (query, donor)
    ]
    renumbered = all_correspondences(
        Chem.RenumberAtoms(query, orders[0]), Chem.RenumberAtoms(donor, orders[1])
    )
    restored = {
        tuple(sorted((orders[0][q], orders[1][d]) for q, d in pairs))
        for pairs in renumbered
    }
    assert restored == direct, "Independent atom renumbering changed valid maps"
    return direct


def compatible(query, donor, mapping):
    return stereo_compatible(
        mapping, stereo_constraints(query), stereo_constraints(donor)
    )


def chemical_examples():
    rng = random.Random(4137)
    # Whole-molecule, chirality-aware RDKit queries retain the stereo that FMCS SMARTS can drop.
    for smi, n_expected in (
        ("C[C@H](F)[C@H](F)C", 1),
        ("C[C@H](F)[C@@H](F)C", 2),
        ("CC/C=C/CC", 2),
        ("CC/C=C\\CC", 2),
    ):
        mol = Chem.MolFromSmiles(smi)
        observed = check_invariance(mol, mol, rng)
        exact = {
            tuple(enumerate(t))
            for t in mol.GetSubstructMatches(mol, useChirality=True, uniquify=False)
        }
        assert len(exact) == n_expected and observed == exact, (smi, observed, exact)
    alkene_e = Chem.MolFromSmiles("CCCCC/C=C/CCCC")
    alkene_z = Chem.MolFromSmiles("CCCCC/C=C\\CCCC")
    partial = check_invariance(alkene_e, alkene_z, rng)
    # The opposite stereoisomer still shares the long connected fragment on one side.
    assert all(len(pairs) == 7 for pairs in partial)
    unspecified = Chem.MolFromSmiles("CCC=CCC")
    known = Chem.MolFromSmiles("CC/C=C/CC")
    assert all(len(pairs) == 6 for pairs in check_invariance(known, unspecified, rng))
    for smi, flipped, center in (
        ("CC[C@H](F)CO", "CC[C@@H](F)CO", 2),
        ("CC[C@](F)(Cl)CO", "CC[C@@](F)(Cl)CO", 2),
    ):
        q, d = Chem.MolFromSmiles(smi), Chem.MolFromSmiles(flipped)
        identity = {i: i for i in range(q.GetNumAtoms())}
        assert not compatible(q, d, identity)
        partial = identity.copy()
        del partial[q.GetAtomWithIdx(center).GetNeighbors()[0].GetIdx()]
        assert compatible(q, d, partial)
        unspecified = Chem.Mol(d)
        Chem.RemoveStereochemistry(unspecified)
        assert compatible(q, unspecified, identity) and compatible(
            unspecified, d, identity
        )

    q = Chem.MolFromSmiles("F/C(Cl)=C(Br)/I")
    transposed = Chem.MolFromSmiles("F/C(Cl)=C(Br)\\I")
    identity = {i: i for i in range(q.GetNumAtoms())}
    assert compatible(q, q, identity) and not compatible(q, transposed, identity)
    stereo_bond = next(b for b in q.GetBonds() if b.GetStereoAtoms())
    start, end = stereo_bond.GetBeginAtomIdx(), stereo_bond.GetEndAtomIdx()
    left, right = stereo_bond.GetStereoAtoms()
    (other_left,) = [
        a.GetIdx()
        for a in q.GetAtomWithIdx(start).GetNeighbors()
        if a.GetIdx() not in (end, left)
    ]
    # Changing one defining substituent reverses the specified alkene relation.
    d = Chem.Mol(q)
    bond = d.GetBondWithIdx(stereo_bond.GetIdx())
    bond.SetStereoAtoms(other_left, right)
    assert not compatible(q, d, identity)
    bond.SetStereo(
        Chem.BondStereo.STEREOZ
        if stereo_bond.GetStereo() == Chem.BondStereo.STEREOE
        else Chem.BondStereo.STEREOE
    )
    assert compatible(q, d, identity)
    partial = {k: v for k, v in identity.items() if k not in (left, other_left)}
    assert compatible(q, transposed, partial)
    nonref = {k: v for k, v in identity.items() if k != left}
    assert compatible(q, d, nonref) and not compatible(q, transposed, nonref)
    Chem.RemoveStereochemistry(d)
    assert compatible(q, d, identity)
    print(
        "Tetrahedral and E/Z controls passed, including atom permutation, exchanged molecules, "
        "changed reference substituent, and partial or unspecified stereo."
    )


def audit_graphs(paths):
    from prepare import ligand_graph, ligand_mapping

    rng = random.Random(29411)
    total = 0
    for panel, folder in zip(("brd4", "t4l", "hcaii"), paths, strict=True):
        root = Path(__file__).resolve().parent / "results" / panel
        records = json.loads((root / "graphs.json").read_text())
        load = ligand_mapping if panel == "hcaii" else ligand_graph
        mols = {lig: load(folder / "prepared", lig)[0] for lig in records}
        for lig, entry in records.items():
            query, donor = mols[lig], mols[entry["donor"]]
            valid = check_invariance(query, donor, rng)
            stored = {tuple(tuple(p) for p in group) for group in entry["maps"]}
            assert stored == valid, (panel, lig, len(stored), len(valid))
            historical = tuple(sorted(tuple(p) for p in entry["single_indices"]))
            assert historical in valid, (
                panel,
                lig,
                "legacy mapping is chemically incompatible",
            )
            print(
                panel,
                lig,
                entry["donor"],
                "stored",
                len(stored),
                "compatible",
                len(valid),
            )
            total += 1
    print(
        f"Verified exchange, independent atom renumbering, and the historical single map for {total} experimental graph pairs."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--panels",
        nargs=3,
        type=Path,
        metavar=("BRD4", "T4L", "HCAII"),
        help="original data directories containing prepared ligand graphs",
    )
    args = parser.parse_args()
    chemical_examples()
    if args.panels:
        audit_graphs(args.panels)
