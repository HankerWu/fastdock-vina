"""Evaluate the repository's full-protein fixed-charge OBC2 comparator on identical top poses."""

import argparse
import itertools
from pathlib import Path

from data import LIGANDS, SEEDS, write
from select_zinc_pose import read_models
from polar_endpoint import PolarEndpoint, atoms


def run(prepared):
    rows = []
    for rec in ("3PTB", "1BTY"):
        protein = atoms(
            (prepared / (rec + ".pdbqt")).read_text().splitlines(), receptor=True
        )
        for name in LIGANDS:
            ligand = atoms((prepared / (name + ".pdbqt")).read_text().splitlines())
            evaluator = PolarEndpoint(protein, ligand)
            for seed, arm in itertools.product(SEEDS, ("native", "shell")):
                path = prepared / f"{rec}_{name}_{seed}_{arm}.pdbqt"
                pose = read_models(path.read_text().splitlines())[1]
                rows.append(
                    dict(
                        receptor=rec,
                        ligand=name,
                        seed=seed,
                        arm=arm,
                        **evaluator.evaluate(atoms(pose)),
                    )
                )
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("prepared", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    write(a.output, run(a.prepared))
