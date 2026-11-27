"""Prepare and dock the three remaining ligands in the original FABP4 competition assay."""

import argparse
import csv
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
from rdkit import Chem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import run_panel
from prepare_inputs import prepare


def read(path):
    with path.open() as h:
        return list(csv.DictReader(h))


MEASURED = read(HERE / "measured.csv")
OLD = ("3P6D", "3P6E")
NEW = ("3P6F", "3P6G", "3P6H")


def prepare_assay(work, original):
    raw = work / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    for pdb in OLD:
        shutil.copytree(original / pdb, work / pdb, dirs_exist_ok=True)
    info = []
    for label in MEASURED[2:]:
        t = time.perf_counter()
        pdb, ccd = label["pdb"], label["ccd"]
        out = work / pdb
        out.mkdir(exist_ok=True)
        pdf = raw / (pdb + ".pdb")
        cif = raw / (ccd + ".cif")
        urls = {
            pdf: f"https://files.rcsb.org/download/{pdb}.pdb",
            cif: f"https://files.rcsb.org/ligands/download/{ccd}.cif",
        }
        hashes = {f.name: prepare.base.fetch(f, url) for f, url in urls.items()}
        lines = pdf.read_text().splitlines()
        refs = prepare.base.dep_records(lines, ccd)
        chem = prepare.base.ligand(cif, refs, out / "ligand.pdbqt")
        mol = Chem.MolFromSmiles(chem["canonical_smiles"])
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if label["id"] in ("L5", "L6"):
            assert [c for _, c in centers] == [{"L5": "R", "L6": "S"}[label["id"]]]
        meta = dict(
            **chem,
            center=np.array(list(refs[0].values())).mean(0).tolist(),
            target="FABP4",
            pdb=pdb,
            receptor_status="query only",
            ca={
                l[22:27].strip() + ":" + l[17:20]: prepare.base.coords(l)
                for l in lines
                if l.startswith("ATOM  ")
                and l[21] == "A"
                and l[16] in " A"
                and l[12:16].strip() == "CA"
            },
        )
        (out / "preparation.json").write_text(json.dumps(meta) + "\n")
        info.append(
            dict(
                id=label["id"],
                source_pdb=urls[pdf],
                source_ccd=urls[cif],
                raw_sha256=hashes,
                ligand_sha256=run_panel.sha(out / "ligand.pdbqt"),
                smiles=chem["canonical_smiles"],
                CIP=centers,
                charge=Chem.GetFormalCharge(mol),
                heavy_atoms=mol.GetNumHeavyAtoms(),
                prepare_seconds=time.perf_counter() - t,
            )
        )
    (work / "sources.json").write_text(json.dumps(info, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "dock"))
    p.add_argument("work", type=Path)
    p.add_argument("--original", type=Path)
    for arm in (*run_panel.ARMS, "counter"):
        p.add_argument("--" + arm, type=Path)
    p.add_argument("--cpu", type=int, default=4)
    p.add_argument("--reverse", action="store_true")
    p.add_argument("--rescore", action="store_true")
    a = p.parse_args()
    if a.action == "prepare":
        prepare_assay(a.work, a.original)
    else:
        run_panel.run(a, systems=[x["pdb"] for x in MEASURED], queries=NEW)
