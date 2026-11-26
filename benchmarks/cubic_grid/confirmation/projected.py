"""Compare relaxed crystal-derived poses to blind top modes using inter plus intra."""
import argparse
import csv
import json
from pathlib import Path
import re
import sys

import numpy as np

TARGETS = {"EPAS1","FKBP1A","RBP4","FABP4","DYRK1A"}


def heavy_atoms(path):
    lines=[r for r in path.read_text().splitlines() if r.startswith(("ATOM  ","HETATM"))
           and r.split()[-1] not in ("H","HD","HS")]
    return lines,np.array([[float(l[i:i+8]) for i in (30,38,46)] for l in lines])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ("fixtures","confirmation","out","binary"):
        ap.add_argument(name,type=Path)
    a=ap.parse_args()
    fixtures,pools,out,binary=[p.resolve() for p in (a.fixtures,a.confirmation,a.out,a.binary)]
    out.mkdir(parents=True,exist_ok=True)
    metadata={}
    for p in sorted(fixtures.glob("*/preparation.json")):
        m=json.loads(p.read_text())
        if m.get("status")=="ready" and m["target"] in TARGETS:
            metadata[m["pdb"]]=m
            link=out/m["pdb"]
            if not link.exists(): link.symlink_to(p.parent,target_is_directory=True)
    (out/"metadata.json").write_text(json.dumps(metadata))
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"paired_crossdocking"))
    from crystal_start import run
    run(out,binary)
    cells=list(csv.DictReader((pools/"confirmation.csv").open()))
    rows=[]
    for s in csv.DictReader((out/"crystal_starts.csv").open()):
        if s["arm"]!="projected" or s["source"]==s["receiver"]: continue
        q,r=s["source"],s["receiver"]
        ligand,lxyz=heavy_atoms(out/"crystal_diagnostic"/f"{q}-{r}-projected.pdbqt")
        rec,rxyz=heavy_atoms(fixtures/r/"receptor.pdbqt")
        dist=np.linalg.norm(lxyz[:,None]-rxyz,axis=2)
        i,j=np.unravel_index(dist.argmin(),dist.shape)
        for c in cells:
            if (q,r)!=(c["query"],c["receiver"]) or c["arm"] not in ("native","cubic"): continue
            first=(pools/c["pool"]).read_text().split("ENDMDL",1)[0]
            blind=float(re.search(r"REMARK INTER \+ INTRA:\s*([-\d.]+)",first)[1])
            rows.append(dict(query=q,receiver=r,seed=c["seed"],arm=c["arm"],blind_first=c["first"],
                             projected_rmsd=s["rmsd"],blind_raw=blind,projected_raw=s["inter_plus_intra"],
                             gap=float(s["inter_plus_intra"])-blind,
                             closest_A=dist[i,j],ligand_atom=ligand[i][6:26].strip(),protein_atom=rec[j][12:26].strip()))
    with (out/"cross_diagnosis.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=rows[0],lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__=="__main__":
    main()
