"""Decompose returned Top-1 pose error after selecting the best receptor-frame chemical map."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from input_shape import fitted_rmsd


def first_pose(path, serial):
    xyz = {}
    for line in path.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM  ", "HETATM")):
            xyz[int(line[6:11])] = [float(line[i:i+8]) for i in (30,38,46)]
    return np.array([xyz[i] for i in serial])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", type=Path)
    p.add_argument("results", type=Path)
    p.add_argument("csv", type=Path)
    opts = p.parse_args()
    rows = []
    for run in csv.DictReader((opts.results/"structures.csv").open()):
        pdb = run["pdb"]
        meta = json.loads((opts.inputs / pdb / "preparation.json").read_text())
        xyz = first_pose(opts.results/"raw"/(run["key"]+".pdbqt"),meta["serial"])
        refs = np.asarray(meta["mapped"])
        d2 = np.mean(np.sum((xyz-refs)**2,axis=2),axis=1)
        best = np.argmin(d2)
        ref = refs[best]
        pose = np.sqrt(d2[best])
        if abs(pose-float(run["top1_A"]))>1e-6:
            raise ValueError("Receptor-frame pose does not reproduce the external evaluation")
        shift = np.linalg.norm(xyz.mean(axis=0)-ref.mean(axis=0))
        internal = fitted_rmsd(xyz,ref)
        rotation = np.sqrt(max(0,d2[best]-shift**2-internal**2))
        rows.append(dict(key=run["key"],pdb=pdb,arm=run["arm"],seed=run["seed"],cap=run["cap"],
                         rmsd_A=f"{pose:.6f}",centroid_A=f"{shift:.6f}",
                         residual_after_rigid_fit_A=f"{internal:.6f}",rotation_component_A=f"{rotation:.6f}"))
    with opts.csv.open("w",newline="") as f:
        out=csv.DictWriter(f,rows[0],lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


if __name__=="__main__":
    main()
