"""Freeze independent Fig. 4 data and use the existing chemical preparation recipe."""
import argparse
import csv
import hashlib
import sys
import tarfile
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "kinase_state_affinity/source"))
from prepare import ligand

HERE = Path(__file__).resolve().parent
PAPER = "10.1021/ja043475p"
# Independent numeric transcription of the original Figure 4 in compound-number order.
FIG = [55, 6.7, 26, 11, 13, 2.5, 3.6, 18, 29, 7.9, 4, 1.1, 17, 7.8,
       10, 7.5, 4.6, 2.7, 4.2, 2.9, .95, 7.3, 3.7, 30, 33, 6.1, 8.2, 3.9,
       1.3, 2.5, .95, 11, 3.8, 1.2, 3.5, 4.6, 3, .48, 4.4, 9.6, .32, 8.9,
       5.7, 4.2, 2.1, 2.3, 1.6, 1.6, 1.6, 9.1, 2.6, None, 2.8, 2.5, 1.9, 3.4]


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def prep(job):
    work, r = job
    start = time.perf_counter()
    ligand(work, {"smiles": r["smiles"], "cid": r["compound"]})
    out = work / "prepared" / (r["compound"] + ".pdbqt")
    return {"compound": r["compound"], "seconds": time.perf_counter()-start,
            "sha256": hashlib.sha256(out.read_bytes()).hexdigest()}


def main(work, jobs):
    work.mkdir(parents=True, exist_ok=True)
    with tarfile.open(ROOT / "galectin_affinity/evidence.tar.xz") as t:
        t.extractall(work, members=[x for x in t if x.name in
            {"raw/GFA.csv", "panel.csv", "prepared/1KJL.box", "prepared/3T1L.box",
             "prepared/1KJL.pdbqt", "prepared/3T1L.pdbqt"}], filter="data")
    with (work / "panel.csv").open() as f:
        old = [Chem.MolFromSmiles(r["smiles"]) for r in csv.DictReader(f)]
    overlap = {Chem.MolToSmiles(m) for m in old}
    old_scaffolds = {MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in old}
    selected = []
    with (work / "raw/GFA.csv").open(encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            if r["Citation"] != PAPER or r["Gal-3 Kd operator"] != "=":
                continue
            idx = int(r["Name in source"])
            kd = float(r["Gal-3 Kd value (μM)"])
            if FIG[idx-4] != kd:
                raise ValueError((idx,kd,FIG[idx-4]))
            if idx in (4, 12, 44, 55):
                continue
            mol = Chem.MolFromSmiles(r["SMILES"])
            graph = Chem.MolToSmiles(mol)
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
            if graph in overlap:
                raise ValueError("overlap " + str(idx))
            selected.append({"compound": str(idx), "gfa_id":r["GFA_ID"], "kd_uM":kd,
                             "heavy":mol.GetNumHeavyAtoms(), "charge":Chem.GetFormalCharge(mol),
                             "smiles":graph,"scaffold":scaffold,"prior_scaffold":scaffold in old_scaffolds})
    selected.sort(key=lambda r:int(r["compound"]))
    assert len(selected)==len({r["smiles"] for r in selected})==47
    write(work/"manifest.csv",selected)
    if jobs:
        with ProcessPoolExecutor(jobs) as pool:
            costs = list(pool.map(prep,((work,r) for r in selected)))
        write(work/"prep.csv",costs)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work",type=Path)
    parser.add_argument("--prepare-jobs",type=int,default=0)
    args=parser.parse_args()
    main(args.work,args.prepare_jobs)
