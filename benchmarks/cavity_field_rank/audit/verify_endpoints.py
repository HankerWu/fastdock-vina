"""Recount primary endpoints from the frozen poses, raw stages, published Ki table and CLI outputs."""
import argparse
import csv
import gzip
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau, spearmanr
from direct_space import deposited_frames, distances, equivalent_targets, table

def read(path):
    with (gzip.open(path,"rt") if path.suffix==".gz" else path.open()) as f:
        return list(csv.DictReader(f))

def recount(work, fk, pose_tables, fk_tables, out):
    data=np.load(out/"direct.npz")
    frames=deposited_frames(work/"audit_frames")
    inputs=json.loads((work/"inputs.json").read_text())["cells"]
    stages={(r["stem"],r["stage"],r["scope"],r["index"]):float(r["rmsd_A"]) for r in read(pose_tables/"stages.csv")}
    nstage=len(stages)
    worst_stage=0
    for stem in sorted({k[0] for k in stages}):
        drug=stem.split("_")[0]
        ids,ref=equivalent_targets(drug,work/inputs[drug]["ligand"],frames)
        lines=(work/"raw"/(stem+".stages")).read_text().splitlines()
        serial=list(map(int,lines[0][1:].split()))
        for line in lines[1:]:
            tag,tr,idx,text=line.split(maxsplit=3)
            key=(stem,tag,tr,idx)
            x=np.fromstring(text,sep=" ").reshape(-1,3)[[serial.index(i) for i in ids]]
            d=distances(x[None],ref)[0]
            worst_stage=max(worst_stage,abs(d-stages.pop(key)))
    assert not stages
    samples=read(pose_tables/"candidates.csv.gz")
    worst_cand=max(abs(data[r["case"]+"_true"][int(r["index"]),0]-float(r["rmsd_A"])) for r in samples)
    modes=read(out/"independent_modes.csv")
    pooled={}
    for r in modes:
        key=r["case"]+"_"+r["arm"]
        group=pooled.setdefault(key,[])
        assert int(r["rank"])==len(group)+1
        group.append(float(r["rmsd_A"]))
    worst_endpoint=0
    for r in read(pose_tables/"cells.csv"):
        g=pooled.pop(r["stem"],[])
        assert len(g)==int(r["n"])
        if g:
            for k,v in zip(("top1_A","best5_A","oracle_A","regret_A"),(g[0],min(g[:5]),min(g),g[0]-min(g)),strict=True):
                worst_endpoint=max(worst_endpoint,abs(float(r[k])-v))
    assert not pooled and max(worst_endpoint,worst_cand,worst_stage)<1e-8

    # Transcribed independently from Table II, p. 9927 (Holt et al. 1993), nM.
    low=np.array([186,110,12,250,10,7,300,300])
    high=low.copy()
    high[6]=600
    definitive=[(i,j) for i,j in itertools.combinations(range(8),2) if low[j]>=3*high[i] or low[i]>=3*high[j]]
    assert len(definitive)==15
    assay={(r["arm"],r["receiver"]):r for r in read(fk_tables/"assay.csv")}
    evidence=[]
    for arm in sorted({k[0] for k in assay}):
        scores=[]
        for receiver in ("1FKG","1FKH"):
            group=[]
            for n in range(4,12):
                lines=(fk/"pools"/arm/f"{receiver}_{n}.pdbqt").read_text().splitlines()
                group.append(float(next(l for l in lines if l.startswith("REMARK VINA RESULT:")).split()[3]))
            scores.append(group)
        for receiver,pred in zip(("1FKG","1FKH","mean"),[*np.array(scores),np.mean(scores,axis=0)],strict=True):
            points=low==high
            rho=spearmanr(low[points],pred[points]).statistic
            tau=kendalltau(low[points],pred[points]).statistic
            credit=0
            for i,j in definitive:
                direction=np.sign((pred[j]-pred[i])*(low[j]-low[i]))
                credit+=(1+direction)/2
            r=assay.pop((arm,receiver))
            np.testing.assert_allclose([credit,rho,tau,pred[6]-pred[4],pred[7]-pred[5]],
                [float(r[k]) for k in ("correct","spearman","kendall_b","cliff8_10","cliff9_11")],atol=1e-10)
            evidence.append(dict(arm=arm,receiver=receiver,correct=credit,denominator=len(definitive),
                                 rho=rho,tau_b=tau,score10_minus8=pred[6]-pred[4],score11_minus9=pred[7]-pred[5]))
    assert not assay
    table(out/"independent_assay.csv",evidence)
    (out/"verification.json").write_text(json.dumps(dict(candidates=len(samples),stages=nstage,modes=len(modes),
        max_candidate_delta_A=worst_cand,max_stage_delta_A=worst_stage,max_endpoint_delta_A=worst_endpoint,
        ki_table_page=9927,source_doi="10.1021/ja00075a008",assay_checks=len(evidence)),indent=2)+"\n")

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("pose_work","fkbp_work","pose_tables","fkbp_tables","out"):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    recount(a.pose_work,a.fkbp_work,a.pose_tables,a.fkbp_tables,a.out)
