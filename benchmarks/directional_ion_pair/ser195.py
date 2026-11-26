"""Test the Ser195 contact proposed by explicit-solvent FEP against frozen trypsin poses."""

import argparse
import itertools
import numpy as np
from pathlib import Path
from rdkit import Chem
from evaluate import pose_models
from contrasts import EXTRA
from experiment import ARMS, SEEDS, write
from score_physical_pool import serial_index_pairs
from make_ion_pair_anchors import xyz
from trypsin import RECEPTORS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    a = parser.parse_args()
    work = a.work / "trypsin"
    results = []
    for rec, lig in itertools.product(RECEPTORS, ("NH2", "OCH3")):
        lines = (work / f"{lig}.pdbqt").read_text().splitlines()
        smiles = next(
            line.split()[2]
            for line in lines
            if line.startswith("REMARK SMILES ")
            and not line.startswith("REMARK SMILES IDX ")
        )
        mol = Chem.MolFromSmiles(smiles)
        ((_, sub),) = mol.GetSubstructMatches(
            Chem.MolFromSmarts("[c]-[N;H2;+0]" if lig == "NH2" else "[c]-[O;X2]")
        )
        serial = dict(
            (i - 1, s) for i, s in serial_index_pairs(lines, "REMARK SMILES IDX ")
        )[sub]
        (ox,) = [
            np.asarray(xyz(line))
            for line in (work / f"{rec}.pdbqt").read_text().splitlines()
            if line.startswith("ATOM")
            and line[17:20] == "SER"
            and line[21] == "A"
            and int(line[22:26]) == 195
            and line[12:16].strip() == "OG"
        ]
        for seed, arm in itertools.product((*SEEDS, *EXTRA), ARMS):
            models = pose_models(work / f"{rec}_{lig}_{seed}_{arm}.pdbqt")
            dist = [float(np.linalg.norm(x[serial] - ox)) for _, x in models]
            results.append(
                dict(
                    receptor=rec,
                    ligand=lig,
                    seed=seed,
                    arm=arm,
                    model1_A=dist[0],
                    best_A=min(dist),
                    within_3p5=sum(v <= 3.5 for v in dist),
                    modes=len(dist),
                )
            )
    write(a.output, results)


if __name__ == "__main__":
    main()
