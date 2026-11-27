"""Post hoc receptor-frame geometry, exact injected-mode lineage and published pools."""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

from pose_runs import COUNT_KEYS, OLD, digest, table

sys.path.insert(0, str(OLD / "mutation_scan"))
import run as scan
import quality_cost


def rows(f):
    with f.open() as stream:
        return list(csv.DictReader(stream))


def rms(x, refs):
    return np.sqrt(np.square(x - refs).sum(axis=-1).mean(axis=-1)).min()


def inspect_fork(p,stem,cell,tr,step):
    f=p/"raw"/(stem+".stages")
    records=f.read_text().splitlines()
    ids=list(map(int,records[0][1:].split()))
    refs=cell["refs"][:,[list(cell["serials"]).index(i) for i in ids]]
    stages=[]
    positions={}
    for line in records[1:]:
        tag,scope,index,*vals=line.split()
        scope,index=int(scope),int(index)
        x=np.array(vals,float).reshape(-1,3)
        positions.setdefault(tag,{})[(scope,index)]=x
        stages.append(dict(stem=stem,stage=tag,scope=scope,index=index,rmsd_A=rms(x,refs)))
    log=(p/"raw"/(stem+".inject")).read_text().splitlines()
    second=next((s for s in log if s.startswith("second ")),None)
    kept=int(second.split()[-1]) if second else 0
    selected=positions["selected"][tr,step]
    returned=positions["hunt"][tr,step]
    offered=positions["offered"]
    final_pool=scan.replay.geom.mode_coords(p/"raw"/(stem+".pdbqt"),ids)
    lineage=dict(inserted_id="",first_eviction_step="",end_trajectory_rank="",
                 pooled_index="",fin_input_index="",direct_release_rank="")
    saved=positions.get("second",{}).get((tr,step))
    history=[]
    recorded_id=None
    was_resident=False
    for line in (p/"raw"/(stem+".history")).read_text().splitlines():
        k,e,cap,threshold,*members=line.split()
        k, members=int(k),list(map(int,members))
        member_errs=[rms(offered[tr,idx],refs) for idx in members]
        if k==step and kept:
            recorded_id=next(idx for idx in members if np.array_equal(offered[tr,idx],saved))
            lineage["inserted_id"]=recorded_id
        resident=recorded_id in members
        if was_resident and not resident and lineage["first_eviction_step"]=="":
            lineage["first_eviction_step"]=k
        was_resident=resident
        history.append(dict(stem=stem,step=k,offered_A=rms(offered[tr,k],refs),score=float(e),
                            capacity=int(cap),dedup_radius_A=float(threshold),
                            resident_ids=" ".join(map(str,members)),
                            resident_oracle_A=min(member_errs,default=""),
                            source_resident=int(resident)))
    if saved is not None and kept:
        def ranks(tag,scope):
            return [i for (sc,i),x in positions.get(tag,{}).items() if sc==scope and np.array_equal(x,saved)]
        for key,tag,scope in (("end_trajectory_rank","end_trajectory",tr),
                              ("pooled_index","pooled",0),("fin_input_index","fin_input",0)):
            found=ranks(tag,scope)
            if found:
                lineage[key]=found[0]
        if lineage["fin_input_index"]!="":
            x=positions["fin_output"][(0,lineage["fin_input_index"])]
            # PDBQT gives coordinates to three decimal places.
            rank=next((i for i,y in enumerate(final_pool,1) if np.max(abs(x-y))<0.000501),"")
            lineage["direct_release_rank"]=rank
    return stages, history, dict(stem=stem,selected_A=rms(selected,refs),hunt_A=rms(returned,refs),
                                  hunt_evals=int(next(s for s in log if s.startswith("hunt_evals")).split()[1]),
                                  accepted=int("accepted" in log),second=int(second is not None),inserted=kept,
                                  **lineage)


def main(p):
    cells=scan.config(p/"aurora",p/"brd",p)
    costs=rows(p/"control_cost.csv")+rows(p/"branch_cost.csv")
    pools, modes, stages, events, forkstats=[],[],[],[],[]
    for cost in costs:
        stem=cost["stem"]
        name,seed,arm=stem.split("_",2)
        group="development" if seed=="61231" else "independent_seed"
        typed_arm=arm
        if arm.startswith("t"):
            t,s,typed_arm=arm.split("_")
            stage,event,chain=inspect_fork(p,stem,cells[name],int(t[1:]),int(s[1:]))
            stages+=stage
            events+=event
            forkstats.append(chain)
        pool=p/"raw"/(stem+".pdbqt")
        assert digest(pool)==cost["pool_sha256"]
        m=scan.metrics(cells[name],pool)
        errors=np.array([r for r,_ in m])
        for rank,(r,score) in enumerate(m,1):
            modes.append(dict(stem=stem,drug=name,seed=seed,arm=typed_arm,rank=rank,rmsd_A=r,score=score))
        n=len(errors)
        pools.append(dict(stem=stem,group=group,drug=name,seed=seed,arm=arm,method=typed_arm,
                          modes=n,top1_A=errors[0] if n else "",best5_A=min(errors[:5],default=""),
                          oracle_A=min(errors,default=""),regret_A=errors[0]-min(errors) if n else "",
                          first_le2_rank=next((j for j,r in enumerate(errors,1) if r<=2),""),
                          first_le3_rank=next((j for j,r in enumerate(errors,1) if r<=3),""),
                          **{k:cost[k] for k in ("wall_s","cpu_s",*COUNT_KEYS)}))
    table(p/"cells.csv",pools)
    quality_cost.summary(p/"aurora",p/"brd",p)
    table(p/"modes.csv",modes)
    table(p/"stages.csv",stages)
    table(p/"residency.csv",events)
    table(p/"forks.csv",forkstats)
    dist=[]
    for group,method in sorted({(r["group"],r["method"]) for r in pools}):
        sub=[r for r in pools if r["group"]==group and r["method"]==method]
        ok=[r for r in sub if r["modes"]]
        top=np.array([r["top1_A"] for r in ok])
        dist.append(dict(group=group,method=method,n=len(sub),empty=len(sub)-len(ok),
                         **{k+"_mean":np.mean([r[k] for r in ok]) for k in ("top1_A","best5_A","oracle_A","regret_A")},
                         **dict(zip(("q25","median","q75","p90","p95"),np.percentile(top,[25,50,75,90,95]),strict=True)),
                         top_gt5=int(sum(top>5)),top_ge8=int(sum(top>=8)),top_gt10=int(sum(top>10)),
                         wall_s=sum(float(r["wall_s"]) for r in sub),cpu_s=sum(float(r["cpu_s"]) for r in sub),
                         actual_evals=sum(int(r[k]) for r in sub for k in ("hunt_evals","second_evals","explicit_evals")),
                         late_prescore=sum(int(r["late_prescore"]) for r in sub)))
    table(p/"distributions.csv",dist)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work",type=lambda v:Path(v).resolve())
    main(parser.parse_args().work)
