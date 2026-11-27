"""Frozen-conformer crystal-placement controls for the six high-budget oracle poses."""
import argparse
import csv
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "macrocycle_proposal" / "external"))
from external_score import poses, write


def transform(xyz, refs):
    center = xyz.mean(axis=0)
    target_center = refs[0].mean(axis=0)
    x0 = xyz - center
    trials = []
    for target in refs:
        y0 = target - target_center
        u, _, vt = np.linalg.svd(x0.T @ y0)
        matrix = (u * [1, 1, np.sign(np.linalg.det(u @ vt))]) @ vt
        residual = np.mean(np.sum((x0 @ matrix - y0)**2, axis=1))
        trials.append((residual, matrix))
    _, matrix = min(trials,key=lambda t:t[0])
    return center, target_center, matrix


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "results", "engine", "output"):
        p.add_argument(name,type=Path)
    a = p.parse_args()
    a.output.mkdir(exist_ok=True, parents=True)
    with (a.results / "models.csv").open() as f:
        models = list(csv.DictReader(f))
    rows = []
    for path in sorted((a.results/"raw").glob("*_fixed.pdbqt")):
        key = path.stem
        pdb = key[:4]
        group = [r for r in models if r["key"] == key]
        best = min(group, key=lambda r:float(r["rmsd_A"]))
        rank = int(best["rank"])
        block = re.findall(r"^MODEL .*?^ENDMDL",path.read_text(),re.M|re.S)[rank-1]
        raw = list(poses(path))[rank-1]
        d = a.inputs / pdb
        meta = json.loads((d/"preparation.json").read_text())
        xyz = np.array([raw[i][1] for i in meta["serial"]])
        refs = np.array(meta["mapped"])
        center, target_center, rot = transform(xyz, refs)
        receptor = next(poses(d/"receptor.pdbqt"))
        tree = cKDTree([x for t,x in receptor.values() if t not in ("H", "HD")])
        for arm, r, target in (("native",np.eye(3),center),
                               ("translate",np.eye(3),target_center),
                               ("rigid",rot,target_center)):
            dest = a.output / f"{key}_{arm}.pdbqt"
            lines = []
            for line in block.splitlines():
                if line.startswith(("ATOM  ","HETATM")):
                    x = (raw[int(line[6:11])][1]-center)@r + target
                    line = line[:30]+ "".join(f"{v:8.3f}" for v in x) + line[54:]
                lines.append(line)
            dest.write_text("\n".join(lines)+"\n")
            log = subprocess.check_output([a.engine.resolve(),"--score_only","--pose_coordinates",dest,
                "--ligand",d/"ligand.pdbqt","--receptor",d/"receptor.pdbqt","--config",a.results/(pdb+".box"),
                "--cpu","1","--seed","458921","--unbound_energy","0"],text=True,stderr=subprocess.STDOUT)
            dest.with_suffix(".log").write_text(log)
            terms = re.findall(r":\s*([-+.\deE]+)\s+\(kcal/mol\)",log)
            assert len(terms)==10, (key,terms)
            pos = (xyz-center)@r+target
            rms = np.sqrt(np.sum((pos-refs)**2,axis=(1,2))/len(xyz)).min()
            nearest = tree.query(pos)[0]
            rows.append(dict(key=key,rank=rank,arm=arm,rmsd_A=rms,nearest_A=nearest.min(),
                             clash_atoms=int((nearest < 2).sum()),energy=float(terms[0]),
                             inter=float(terms[1]),intra=float(terms[4])))
    write(a.output/"placement.csv", rows)
