"""Labelled single-event interventions after all blind predictions have been frozen."""
import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from direct_space import BENCH, deposited_frames, equivalent_targets, distances, table

sys.path.insert(0,str(BENCH/"cavity_triplet"))
from pose_runs import invoke
from pose_metrics import inspect_fork
import json

def forks(pose, blind_results, out, observer):
    (out/"raw").mkdir(parents=True,exist_ok=True)
    frames=deposited_frames(out)
    systems=json.loads((pose/"inputs.json").read_text())["cells"]
    costs, summaries=[],[]
    eligible=defaultdict(list)
    with gzip.open(blind_results/"candidates.csv.gz","rt") as f:
        for row in csv.DictReader(f):
            if float(row["rmsd_A"])<=2:
                eligible[row["case"]].append(row)
    for case, pool in sorted(eligible.items()):
        drug,seed=case.split("_")
        cfg=systems[drug].copy()
        for k in ("receptor","ligand","bonds"):
            cfg[k]=pose/cfg[k]
        ids,refs=equivalent_targets(drug,cfg["ligand"],frames)
        choices=defaultdict(list)
        bank=(pose/"banks"/(case+".bank")).read_text().splitlines()
        for arm in ("shape","true","perm"):
            first=min(pool,key=lambda row:int(row[arm+"_rank"]))
            choices[int(first["index"])].append((arm,first[arm+"_rank"]))
        for i,origins in choices.items():
            stem=f"{case}_{i}"
            root=out/"raw"/stem
            injection=root.with_suffix(".input")
            injection.write_text("1 1 "+bank[i]+"\n")
            env=dict(FD_RIGID_INTERVENTION=injection,FD_RIGID_LOG=root.with_suffix(".inject"),
                     FD_RIGID_COORDS=root.with_suffix(".stages"),FD_RIGID_HISTORY=root.with_suffix(".history"))
            output=invoke(out,costs,stem,observer,cfg,int(seed),env=env)
            stages,_,chain=inspect_fork(out,stem,dict(serials=ids,refs=refs),1,1)
            errors=[]
            for block in output.read_text().split("ENDMDL")[:-1]:
                x={int(l[6:11]):[float(l[j:j+8]) for j in (30,38,46)] for l in block.splitlines() if l.startswith(("ATOM  ","HETATM"))}
                errors.append(distances(np.array([[x[s] for s in ids]]),refs)[0])
            summaries.append(dict(case=case,index=i,rank_in_arms=" ".join(a+"="+r for a,r in origins),
                second_A=next((r["rmsd_A"] for r in stages if r["stage"]=="second"),""),
                top1_A=errors[0],best5_A=min(errors[:5]),oracle_A=min(errors),**chain))
    table(out/"near_forks.csv",summaries)
    table(out/"cost.csv",costs)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for v in ("pose","results","out","observer"):
        p.add_argument(v,type=lambda s:Path(s).resolve())
    a=p.parse_args()
    forks(a.pose,a.results,a.out,a.observer)
