"""Prepare fixed crystal structures and ligand conformers independently of affinity."""

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

os.environ["OPENMM_CPU_THREADS"] = "1"
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from openmm import Platform
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from freeze import RECEPTORS


def receptor(work, pdb, component):
    lines = (work / f"raw/{pdb}.pdb").read_text().splitlines()
    prepared = work / "prepared"
    protein = prepared / f"{pdb}-protein.pdb"
    # Keep native polymer heavy atoms; 2SRC tail pTyr527 is far from the box
    # and is represented by unphosphorylated Tyr in this rigid protein approximation.
    atoms = []
    for line in lines:
        if line[21:22] != "A" or line[16:17] not in (" ", "A"):
            continue
        if (
            line.startswith("HETATM")
            and line[17:20] == "PTR"
            and line[12:16].strip() not in ("P", "O1P", "O2P", "O3P")
        ):
            line = "ATOM  " + line[6:17] + "TYR" + line[20:]
        if line.startswith("ATOM  ") and line[76:78].strip() not in ("H", "D"):
            atoms.append(line[:16] + " " + line[17:])
    protein.write_text("\n".join(atoms) + "\nEND\n")
    fixer = PDBFixer(filename=str(protein), platform=Platform.getPlatformByName("CPU"))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    missing = {
        f"{r.name}:{r.id}": [a.name for a in atoms]
        for r, atoms in fixer.missingAtoms.items()
    }
    fixer.addMissingAtoms(seed=2026)
    with protein.open("w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    with (prepared / f"{pdb}.log").open("w") as log:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "meeko.cli.mk_prepare_receptor",
                "--read_pdb",
                protein,
                "-p",
                prepared / f"{pdb}.pdbqt",
            ],
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    ligand = [
        line
        for line in lines
        if line.startswith("HETATM")
        and line[21] == "A"
        and line[17:20].strip() == component
        and line[16] in " A"
        and line[76:78].strip() not in ("H", "D")
    ]
    assert len(ligand) > 12
    xyz = np.array([[float(l[k : k + 8]) for k in (30, 38, 46)] for l in ligand])
    (prepared / f"{pdb}.box").write_text(
        "".join(
            f"center_{a} = {v:.5f}\nsize_{a} = 26\n"
            for a, v in zip("xyz", xyz.mean(axis=0))
        )
    )
    (prepared / f"{pdb}-missing.json").write_text(json.dumps(missing, indent=2) + "\n")


def ligand(work, row):
    mol = Chem.MolFromSmiles(row["smiles"])
    parts = Chem.GetMolFrags(mol, asMols=True)
    if len(parts) > 1:
        if (
            row["cid"] != "11984591"
            or sorted(Chem.MolToSmiles(x) for x in parts)[-1] != "O=S(=O)(O)c1ccccc1"
        ):
            raise ValueError(f"Unexpected salt components: {row}")
        mol = max(parts, key=lambda m: m.GetNumHeavyAtoms())
    mol = Chem.AddHs(mol)
    opt = AllChem.ETKDGv3()
    opt.randomSeed = 2026
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=8, params=opt))
    optima = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=1500)
    viable = [(energy, i) for i, (status, energy) in zip(ids, optima) if status == 0]
    if not viable:
        raise ValueError(row["cid"])
    _, i = min(viable)
    conf = Chem.Conformer(mol.GetConformer(i))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    (setup,) = MoleculePreparation().prepare(mol)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(err)
    (work / "prepared" / f"{row['cid']}.pdbqt").write_text(pdbqt)


if __name__ == "__main__":
    work = Path(sys.argv[1])
    work.joinpath("prepared").mkdir(exist_ok=True)
    with (work / "panel.csv").open() as f:
        rows = {r["cid"]: r for r in csv.DictReader(f)}
    mode = sys.argv[2]
    if mode == "ligand":
        for row in rows.values():
            ligand(work, row)
    else:
        for states in RECEPTORS.values():
            for pdb, component in states:
                if mode in ("receptor", pdb):
                    receptor(work, pdb, component)
