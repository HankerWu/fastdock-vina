"""Score raw fixed pose pools, then evaluate against deposited receptor-frame poses."""

import argparse
import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolAlign

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

from rank_halogen_pool import carbonyls, corrections, models, rank_pool


def write(path, rows):
    with path.open("w") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def analyze(work, output):
    output.mkdir(exist_ok=True, parents=True)
    poses, cells, crystals = [], [], []
    for folder in sorted(work.iterdir()):
        if not (folder / "ready").exists() or not list(folder.glob("poses-*.pdbqt")):
            continue
        ref = Chem.SDMolSupplier(str(folder / "crystal.sdf"), removeHs=True)[0]
        receptor = folder / "receptor.pdbqt"
        radial, directed = corrections(
            ref, ref.GetConformer().GetPositions(), carbonyls(receptor)
        )
        crystals.append(
            dict(case=folder.name, radial_delta=radial, directional_delta=directed)
        )
        for path in sorted(folder.glob("poses-*.pdbqt")):
            rows = rank_pool(receptor, path)
            seed = int(path.stem.split("-")[-1])
            for row, (model, score, mol, xyz) in zip(rows, models(path), strict=True):
                if model != row["model"] or score != row["native"]:
                    raise ValueError("inconsistent model order")
                c = Chem.Conformer(mol.GetNumAtoms())
                c.SetPositions(xyz)
                mol.AddConformer(c)
                row["rmsd"] = rdMolAlign.CalcRMS(mol, ref, maxMatches=100000)
                row["case"] = folder.name
                row["seed"] = seed
            poses.extend(rows)
            cell = dict(
                case=folder.name,
                seed=seed,
                models=len(rows),
                oracle_rmsd=min(r["rmsd"] for r in rows),
            )
            for arm in ("native", "radial", "directional"):
                (selected,) = [r for r in rows if r[arm + "_rank"] == 1]
                cell[arm + "_model"] = selected["model"]
                cell[arm + "_rmsd"] = selected["rmsd"]
            cells.append(cell)
    write(output / "poses.csv", poses)
    write(output / "cells.csv", cells)
    write(output / "crystals.csv", crystals)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    analyze(args.work, args.output)
