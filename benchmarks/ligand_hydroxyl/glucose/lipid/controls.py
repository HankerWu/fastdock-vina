"""Post-discovery controls for the distal lipid's non-CCD stereochemistry and score terms."""
import argparse
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from restore import dock, write_csv
from paired import atoms, cli_score, crystal_block


def controls(w,o,binary):
    t=o/"near_only"
    t.mkdir(exist_ok=True)
    lines=(o/"4ZW9-OLC.pdbqt").read_text().splitlines()
    keep=[l for l in lines if not (l[17:20]=="OLC" and int(l[22:26])==503)]
    if len(keep)+27!=len(lines):
        raise ValueError("Failed to remove exactly one OLC")
    (t/"4ZW9-OLC.pdbqt").write_text("\n".join(keep)+"\n")
    dock(w,t,binary)
    meta=json.loads((w/"metadata.json").read_text())
    cells=[]
    for path in sorted((t/"runs").glob("*.pdbqt")):
        case,cap,seed,_=path.stem.split("_")
        m=meta[case]
        rows=[]
        for b in path.read_text().split("ENDMDL"):
            if b.strip():
                v=atoms(b)
                d=np.asarray(m["mapped"])-np.array([v[i].xyz for i in m["serial"]])
                rows.append(float(np.sqrt((d*d).sum(axis=2).mean(axis=1)).min()))
        cells.append(dict(case=case,cap=cap,seed=seed,modes=len(rows),top1=rows[0],top5=min(rows[:5]),oracle=min(rows),
                          identical_full_output=(o/"runs"/path.name).read_bytes()==path.read_bytes()))
    write_csv(o/"distant_control.csv",cells)
    # Native weight ablation at both deposited poses and each blind native rank 1.
    scores=[]
    cases=list(csv.DictReader((o/"poses.csv").open()))
    with TemporaryDirectory() as tmp:
        pose=Path(tmp)/"one.pdbqt"
        for row in cases:
            if row["rank"]!="1":
                continue
            code,cap,seed,pool=[row[k] for k in ("code","cap","seed","pool")]
            lig=w/"4ZW9"/(code+".pdbqt")
            if pool=="crystal":
                m=meta["4ZW9-"+code]
                b=crystal_block(lig.read_text(),m["mapped"][0],m["serial"])
            else:
                suff="pro" if pool=="dry" else "olc"
                p=(w if pool=="dry" else o)/"runs"/f"4ZW9-{code}_{cap}_{seed}_{suff}.pdbqt"
                b=p.read_text().split("ENDMDL",1)[0]+"ENDMDL\n"
            pose.write_text(b)
            fields={}
            for label,flag in [("no_hb","--weight_hydrogen"),("no_rep","--weight_repulsion")]:
                dry=cli_score(binary,w/"4ZW9/receptor.pdbqt",lig,pose,[flag,"0"])
                olc=cli_score(binary,o/"4ZW9-OLC.pdbqt",lig,pose,[flag,"0"])
                fields[label+"_raw_delta"]=olc[0]-dry[0]
            near=cli_score(binary,t/"4ZW9-OLC.pdbqt",lig,pose)
            scores.append(dict(code=code,pool=pool,cap=cap,seed=seed,rank=row["rank"],
                          near_raw_delta=near[0]-float(row["dry_inter"]),full_raw_delta=float(row["raw_inter_delta"]),**fields))
    write_csv(o/"weight_controls.csv",scores)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("sugar",type=Path)
    p.add_argument("out",type=Path)
    p.add_argument("binary",type=Path)
    a=p.parse_args()
    controls(a.sugar.resolve(),a.out.resolve(),a.binary.resolve())
