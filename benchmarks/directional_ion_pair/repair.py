"""Check a strained fixed input against independent free ligand conformers before pairing."""

import argparse
import itertools
from pathlib import Path

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

from experiment import ARMS, ROOT, SEEDS, prepare, run, table, write

LIG = "CHEMBL369480"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("out", "fjs", "p16", "binary"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    row = next(
        r
        for r in table(ROOT / "benchmarks/affinity_partition/fxa_series.csv")
        if r["ligand"] == LIG
    )
    smiles = row["docked_smiles"]
    attempts = []
    best = None
    for seed in range(2026, 2042):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        par = AllChem.ETKDGv3()
        par.randomSeed = seed
        if AllChem.EmbedMolecule(mol, par):
            continue
        status = AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        props = AllChem.MMFFGetMoleculeProperties(mol)
        score = AllChem.MMFFGetMoleculeForceField(mol, props).CalcEnergy()
        coords = mol.GetConformer().GetPositions()
        d = max(
            np.linalg.norm(coords[b.GetBeginAtomIdx()] - coords[b.GetEndAtomIdx()])
            for b in mol.GetBonds()
            if b.GetBeginAtom().GetAtomicNum() > 1 and b.GetEndAtom().GetAtomicNum() > 1
        )
        attempts.append(
            dict(
                seed=seed,
                converged=int(status == 0),
                mmff_kcal=score,
                longest_heavy_bond=d,
            )
        )
        if status == 0 and (best is None or score < best[0]):
            best = score, mol
    write(a.out / "conformer_selection.csv", attempts)
    (setup,) = MoleculePreparation().prepare(best[1])
    text, ok, why = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(why)
    ligand = a.out / f"{LIG}.pdbqt"
    ligand.write_text(text)
    inputs = {}
    for rec, directory in (("1FJS", a.fjs), ("2P16", a.p16)):
        prepare(a.out, rec, directory / f"{rec}.pdb", "A", [ligand])
        inputs[rec, LIG] = ligand, directory / f"{rec}.pdbqt", directory / "box.txt"
    run(
        a.binary.resolve(),
        a.out,
        inputs,
        itertools.product(("1FJS", "2P16"), (LIG,), (*SEEDS, 161803, 141421), ARMS),
        2,
    )


if __name__ == "__main__":
    main()
