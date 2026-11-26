"""Audit bound coordinate RMSDs, nitrogen handedness and graph correspondence invariance."""

import argparse
import csv
import json
import random
import sys
from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory

from prepare import (
    ROOT,
    all_correspondences,
    declared_ca_graph,
    ligand_graph,
    ligand_mapping,
)
from rdkit import Chem

sys.path.insert(0, str(ROOT / "benchmarks/cutoff_switch"))
from verify_cutoff_poses import verify as verify_rmsd
from verify_hcaii import verify as verify_signed


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def save(path, rows):
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def audit(panel, data, runs):
    graphs = json.loads((runs / "graphs.json").read_text())
    load = ligand_mapping if panel == "hcaii" else ligand_graph
    structures = {lig: load(data / "prepared", lig)[0] for lig in graphs}
    if panel == "hcaii" and all(
        v.get("nitrogen") == "declared" for v in graphs.values()
    ):
        structures = {
            lig: declared_ca_graph(data / "prepared", lig, mol)
            for lig, mol in structures.items()
        }
    rng = random.Random(20261126)
    for lig, entry in graphs.items():
        query, donor = structures[lig], structures[entry["donor"]]
        orders = [
            rng.sample(range(m.GetNumAtoms()), m.GetNumAtoms()) for m in (query, donor)
        ]
        maps = all_correspondences(
            *(
                Chem.RenumberAtoms(m, order)
                for m, order in zip((query, donor), orders, strict=True)
            )
        )
        recovered = {
            tuple(sorted((orders[0][a], orders[1][b]) for a, b in mapping))
            for mapping in maps
        }
        expected = {tuple(tuple(p) for p in mapping) for mapping in entry["maps"]}
        if recovered != expected:
            raise ValueError(
                f"Graph embeddings change under atom permutation for {lig}"
            )
    print(
        f"{len(structures)} query/donor graph ensembles invariant to independent atom permutations"
    )
    modes = defaultdict(list)
    for row in table(runs / "modes.csv"):
        modes[tuple(row[k] for k in ("receptor", "ligand", "seed", "arm"))].append(row)
    cells = []
    with TemporaryDirectory() as folder:
        out = Path(folder)
        (out / "raw").symlink_to(data.resolve() / "raw")
        (out / "runs").mkdir()
        for (rec, lig, seed, arm), rows in modes.items():
            name = f"{rec}_{lig}_{seed}_{arm}.pdbqt"
            src = runs.resolve() / name
            if panel == "hcaii":
                (out / name).symlink_to(src)
                for row in rows:
                    cells.append(
                        dict(
                            rec=rec,
                            lig=lig,
                            seed=seed,
                            arm=arm,
                            model=row["model"],
                            rmsd=row["rmsd"],
                            sign_ok=row["sign_ok"],
                        )
                    )
            else:
                (out / "runs" / f"{lig}_{rec}_{seed}_{arm}.pdbqt").symlink_to(src)
                cells.append(
                    dict(
                        ligand=lig,
                        receptor=rec,
                        seed=seed,
                        arm=arm,
                        pool_rmsds=",".join(r["rmsd"] for r in rows),
                    )
                )
        if panel == "hcaii":
            save(out / "modes.csv", cells)
            verify_signed(data, out)
        else:
            save(out / "cells.csv", cells)
            verify_rmsd(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("panel", choices=("brd4", "t4l", "hcaii"))
    parser.add_argument("data", type=Path)
    parser.add_argument("runs", type=Path)
    args = parser.parse_args()
    audit(args.panel, args.data, args.runs)
