"""Independently validate receptor-frame cognate RMSDs with PoseBusters SDFs."""

import argparse
import csv
import hashlib
from pathlib import Path

from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdMolAlign


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("dataset", type=Path, help="posebusters_benchmark_set directory")
    p.add_argument("runs", type=Path, help="pose_run.py output")
    p.add_argument("output", type=Path)
    args = p.parse_args()
    with (args.runs / "modes.csv").open() as stream:
        rows = [r for r in csv.DictReader(stream) if r["kind"] == "redock"]
    lookup = {
        (r["query"], r["seed"], r["arm"], int(r["rank"])): float(r["rmsd_A"])
        for r in rows
    }
    reports = []
    for pdb, seed, arm in sorted({k[:3] for k in lookup}):
        (querydir,) = args.dataset.glob(pdb + "_*")
        (sdf,) = querydir.glob("*_ligand.sdf")
        ref = Chem.SDMolSupplier(str(sdf), removeHs=True)[0]
        pose = args.runs / "poses" / f"{pdb}-{pdb}-{seed}-{arm}.pdbqt"
        (mol,) = RDKitMolCreate.from_pdbqt_mol(
            PDBQTMolecule.from_file(str(pose), skip_typing=True)
        )
        mol = Chem.RemoveHs(mol)
        same = Chem.MolToSmiles(mol) == Chem.MolToSmiles(ref)
        for j in range(mol.GetNumConformers()):
            observed = lookup[pdb, seed, arm, j + 1]
            row = {
                "pdb": pdb,
                "seed": seed,
                "arm": arm,
                "rank": j + 1,
                "sdf_sha256": hashlib.sha256(sdf.read_bytes()).hexdigest(),
                "graph_identical": same,
                "pdb_rmsd_A": observed,
                "sdf_rmsd_A": "",
                "abs_delta_A": "",
            }
            if same:
                rmsd = rdMolAlign.CalcRMS(mol, ref, prbId=j, refId=0)
                row.update(sdf_rmsd_A=rmsd, abs_delta_A=abs(rmsd - observed))
            reports.append(row)
    if len(reports) != len(rows):
        raise RuntimeError("incomplete independent audit")
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, list(reports[0]))
        writer.writeheader()
        writer.writerows(reports)


if __name__ == "__main__":
    main()
