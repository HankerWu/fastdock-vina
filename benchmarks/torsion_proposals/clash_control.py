"""Check fixed-pool intraligand geometry in crystal, prepared and returned poses."""

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from meeko import PDBQTMolecule, RDKitMolCreate
from posebusters.modules.distance_geometry import check_geometry
from rdkit import Chem


def read_mol(path):
    mols = RDKitMolCreate.from_pdbqt_mol(
        PDBQTMolecule.from_file(str(path), skip_typing=True)
    )
    if len(mols) != 1 or mols[0] is None:
        raise ValueError(path)
    return mols[0]


def check(mol):
    report = check_geometry(
        mol,
        threshold_clash=0.3,
        threshold_bad_bond_length=0.25,
        threshold_bad_angle=0.25,
    )
    table = report["details"]["clash"]
    entry = table.loc[(table.distance / table.lower_bound).idxmin()]
    i, j = map(int, entry.atom_pair)
    path = Chem.GetShortestPath(mol, i, j)
    return dict(
        a=i,
        b=j,
        atom_pair=entry.atom_types,
        bond_steps=len(path) - 1,
        distance=float(entry.distance),
        bound=float(entry.lower_bound),
        relative=float(entry.distance / entry.lower_bound),
        passes=bool(report["results"]["no_internal_clash"]),
        n_clash=int(report["results"]["number_clashes"]),
    )


def score_pool(item):
    path, rows = item
    mol = read_mol(path)
    if len(rows) != mol.GetNumConformers():
        raise ValueError(path)
    results = []
    for row, conf in zip(rows, mol.GetConformers()):
        model = Chem.Mol(mol)
        model.RemoveAllConformers()
        model.AddConformer(Chem.Conformer(conf), assignId=True)
        results.append(
            {
                k: row[k]
                for k in ("cap", "ligand", "receiver", "seed", "arm", "rank", "rmsd")
            }
            | check(model)
        )
    return results


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def main(work, jobs):
    import json

    files = []
    for c in json.loads((work / "metadata.json").read_text()):
        files.append(work / c / "ligand.pdbqt")
        files.extend(
            work / "crystal_diagnostic" / f"{c}-{c}-{s}.pdbqt"
            for s in ("crystal-input", "projected-input", "crystal", "projected")
        )
    write(
        work / "clash_controls.csv",
        [dict(path=str(p.relative_to(work)), **check(read_mol(p))) for p in files],
    )
    inputs = {}
    for csvfile in sorted(work.glob("modes-*.csv")):
        for r in csv.DictReader(csvfile.open()):
            if r["cap"] not in ("12000", "250000") or r["arm"] == "balanced":
                continue
            key = (
                work
                / f"runs_{r['cap']}"
                / (
                    "_".join(r[k] for k in ("ligand", "receiver", "seed", "arm"))
                    + ".pdbqt"
                )
            )
            inputs.setdefault(key, []).append(r)
    with ProcessPoolExecutor(jobs) as ex:
        groups = list(ex.map(score_pool, inputs.items()))
    rows = [r for g in groups for r in g]
    write(work / "clash_modes.csv", rows)
    out = []
    for g in groups:
        top = g[0]
        passed = next((r for r in g if r["passes"]), None)
        out.append(
            {k: top[k] for k in ("cap", "ligand", "receiver", "seed", "arm")}
            | dict(
                top1_rmsd=top["rmsd"],
                top1_pass=top["passes"],
                selected_rank=passed["rank"] if passed else None,
                selected_rmsd=passed["rmsd"] if passed else None,
                abstained=passed is None,
            )
        )
    write(work / "clash_selection.csv", out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--jobs", type=int, default=4)
    a = p.parse_args()
    main(a.work, a.jobs)
