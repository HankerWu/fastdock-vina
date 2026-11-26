"""Stress-test the availability of a free aryl amidinium with independent ligand conformers."""

import argparse
import itertools
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from experiment import ROOT, table, write
from make_ion_pair_anchors import AMIDINE, ACID


def evaluate(args):
    ligand, smiles, seed = args
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    opt = AllChem.ETKDGv3()
    opt.randomSeed = seed
    if AllChem.EmbedMolecule(m, opt) != 0:
        return []
    _, _, n1, n2 = m.GetSubstructMatch(AMIDINE)
    _, o1, o2 = m.GetSubstructMatch(ACID)
    outputs = []
    for eps in (4, 78.4):
        mol = Chem.Mol(m)
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        props.SetMMFFDielectricConstant(eps)
        force = AllChem.MMFFGetMoleculeForceField(mol, props)
        converged = force.Minimize(maxIts=1000) == 0
        xyz = mol.GetConformer().GetPositions()
        n = xyz[[n1, n2]]
        best = min(
            (np.linalg.norm(n - xyz[list(o)], axis=1) for o in ((o1, o2), (o2, o1))),
            key=lambda v: max(v),
        )
        outputs.append(
            dict(
                ligand=ligand,
                seed=seed,
                dielectric=eps,
                converged=int(converged),
                mmff=force.CalcEnergy(),
                n_o_short=min(best),
                n_o_long=max(best),
                bidentate=int(max(best) <= 3.5),
            )
        )
    return outputs


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--jobs", type=int, default=2)
    a = p.parse_args()
    rows = [
        r
        for r in table(ROOT / "benchmarks/affinity_partition/fxa_series.csv")
        if int(r["included"]) and r["charge"] == "1"
    ]
    args = [
        (r["ligand"], r["docked_smiles"], seed)
        for r, seed in itertools.product(rows, range(2000, 2032))
    ]
    with ProcessPoolExecutor(a.jobs) as pool:
        attempts = [r for result in pool.map(evaluate, args) for r in result]
    write(a.out, attempts)
