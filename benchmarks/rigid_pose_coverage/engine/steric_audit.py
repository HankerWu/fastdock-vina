"""Post hoc, fixed libraries: count ligand–receptor heavy-atom pairs closer than 2 Å."""
import argparse
import csv
import gzip
import json
import tarfile
from pathlib import Path

import numpy as np

from coverage_gate import ROOT, _geometry
from pocket_triangle import Pocket


def audit(parent, new, out):
    exp = ROOT/"benchmarks/rigid_pose_coverage"
    prior = exp/"audit/results"
    near = {}
    with (prior/"near_neighbors.csv").open() as f:
        for r in csv.DictReader(f):
            if int(r["selected"]):
                near.setdefault(r["case"]+"_"+r["arm"], {})[int(r["near_id"])] = float(r["crystal_A"])
    with gzip.open(exp/"pose/results/candidates.csv.gz", "rt") as f:
        for r in csv.DictReader(f):
            if float(r["rmsd_A"]) <= 2:
                near.setdefault(r["stem"], {})[int(r["triangle_id"])] = float(r["rmsd_A"])
    events, witnesses = [], []
    with np.load(prior/"selected_ids.npz") as archived:
        for panel, work in (("historical",parent), ("prospective",new)):
            meta=json.loads((work/"inputs.json").read_text())["cells"]
            suffix = "whole" if panel=="historical" else "coverage"
            for source in sorted((work/"dump").glob("*_"+suffix+"/source.xyz")):
                case=source.parent.name.removesuffix("_"+suffix)
                cfg=meta[case.rsplit("_",1)[0]]
                header, line=source.read_text().splitlines()
                s=np.array(header.split()[1:], int)
                x=np.array(line.split()[3:],float).reshape(-1,3)[np.argsort(s)]
                pocket=Pocket.from_pdbqt(work/cfg["receptor"],cfg["center"],cfg["size"])
                _, anchor, _ = _geometry(pocket,x,np.sort(s))
                for arm in ("edge","feasible_edge","coverage","uniform"):
                    stem=case+"_"+arm
                    if panel=="historical":
                        ids=archived[stem]
                        centers, rots=pocket.fit_triangles(x,x[anchor],pocket.triangles[ids])
                    else:
                        with np.load(work/"dump"/stem/"bank.npz") as b:
                            ids, centers, rots=b["triangle_id"], b["center"], b["matched_rotation"]
                    poses=(x-x.mean(0)) @ rots.transpose(0,2,1) + centers[:,None,:]
                    separation=pocket.tree.query(poses.reshape(-1,3))[0].reshape(len(ids),-1).min(axis=1)
                    hits=near.get(stem,{})
                    member=np.isin(ids,list(hits))
                    assert member.sum()==len(hits)
                    for i in np.flatnonzero(member):
                        witnesses.append(dict(panel=panel,case=case,arm=arm,triangle_id=ids[i],
                                              rmsd_A=hits[ids[i]],min_separation_A=separation[i]))
                    events.append(dict(panel=panel,case=case,arm=arm,all=len(ids),all_lt2=int(sum(separation<2)),
                                       near=len(hits),near_lt2=int(sum((separation<2)&member))))
    # Existing crystal-informed fork: did native refinement resolve the recovered contact?
    stages=[]
    case="3FMZ-A_93241"
    cfg=json.loads((parent/"inputs.json").read_text())["cells"]["3FMZ-A"]
    pocket=Pocket.from_pdbqt(parent/cfg["receptor"],cfg["center"],cfg["size"])
    with tarfile.open(prior/"capture_raw.tar.xz") as tar:
        lines=tar.extractfile("raw/"+case+"_coverage_labelled.stages").read().decode().splitlines()
    for stage in ("selected","hunt","second"):
        words=next(x for l in lines if (x:=l.split())[:3]==[stage,"1","1"])
        coords=np.array(words[3:],float).reshape(-1,3)
        stages.append(dict(case=case,stage=stage,min_separation_A=pocket.tree.query(coords)[0].min()))
    for name, rows in (("steric_counts",events),("steric_near",witnesses),("steric_capture",stages)):
        with (out/(name+".csv")).open("w",newline="") as f:
            w=csv.DictWriter(f,rows[0].keys(),lineterminator="\n")
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("parent_pose","new_pose","output"):
        parser.add_argument(key,type=Path)
    args=parser.parse_args()
    audit(args.parent_pose,args.new_pose,args.output)
