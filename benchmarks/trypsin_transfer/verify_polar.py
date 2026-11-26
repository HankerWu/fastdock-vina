"""Recalculate matched polar-transfer gaps with double-precision OBC2 and direct Coulomb."""

import argparse
import sys
from pathlib import Path

from scipy.spatial.distance import cdist

from data import ROOT, read, write
from select_zinc_pose import read_models

sys.path.insert(0, str(ROOT / "benchmarks/polar_endpoint"))
from independent import gb_energy
from polar_endpoint import atoms


def calculate(prepared, results):
    archived = read(results / "classical.csv")
    rows = []
    for rec in ("1BTY", "3PTB"):
        receptor = atoms(
            (prepared / (rec + ".pdbqt")).read_text().splitlines(), receptor=True
        )
        apo = gb_energy([receptor])
        for ligand in ("NH2", "OCH3"):
            for seed in (271828, 314159, 141421, 161803):
                path = prepared / f"{rec}_{ligand}_{seed}_shell.pdbqt"
                l = atoms(read_models(path.read_text().splitlines())[1])
                direct = (
                    332.06371
                    / 4
                    * (
                        receptor.charges[:, None]
                        * l.charges
                        / cdist(receptor.xyz, l.xyz)
                    ).sum()
                )
                delta = gb_energy([receptor, l]) - apo - gb_energy([l])
                exact = direct + delta
                stored = next(
                    float(r["polar_transfer"])
                    for r in archived
                    if (r["receptor"], r["ligand"], r["seed"], r["arm"])
                    == (rec, ligand, str(seed), "shell")
                )
                rows.append(
                    dict(
                        receptor=rec,
                        ligand=ligand,
                        seed=seed,
                        openmm=stored,
                        float64=exact,
                        absolute_error=abs(exact - stored),
                    )
                )
    write(results / "independent_polar.csv", rows)
    worst = max(r["absolute_error"] for r in rows)
    print(
        len(rows),
        "full-protein double-precision OBC2 checks; max absolute error",
        worst,
    )
    assert worst < 0.025


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("prepared", type=Path)
    p.add_argument("results", type=Path)
    a = p.parse_args()
    calculate(a.prepared, a.results)
