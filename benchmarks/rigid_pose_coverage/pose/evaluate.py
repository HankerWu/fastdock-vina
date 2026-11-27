"""Rebuild every released/selected RMSD from coordinates, CCD chemistry and receiver C-alpha fitting."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import run_panel as run

load = run.parent.importlib.util
spec = load.spec_from_file_location("fragment_eval",run.BENCH/"fragment_anchor/pose/evaluate.py")
historic=load.module_from_spec(spec)
spec.loader.exec_module(historic)
g=historic.old
g.panel.CASES=run.CASES
g.COMPONENTS={"3ANR":"HRM","7AKL":"RK5","3FMZ":"2T1","1RBP":"RTL"}


def evaluate(p, out):
    out.mkdir(parents=True,exist_ok=True)
    refs=g.reference_inputs(p)
    costs={r["stem"]:r for f in ("run_cost","control_cost","ladder","ladder_counters") for r in run.old.csvrows(p/(f+".csv"))}
    matched={r["case"]:r for r in run.old.csvrows(p/"cost_selected.csv")}
    cells,modes,quality,branches,stages,history,dist,paired=[],[],[],[],[],[],[],[]
    native=("pro","pro-high","official","official-high","pro-cost")
    with gzip.open(out/"candidates.csv.gz","wt",newline="") as stream:
        cw=csv.writer(stream)
        cw.writerow(("stem","index","triangle_id","edge_residual_A2","rmsd_A","field_score","field_rank"))
        for name,(query,rec,seeds) in run.CASES.items():
            cell,mol,src,tree=refs[name]
            for seed in seeds:
                base=f"{name}_{seed}"
                for arm in (*run.ARMS,*native):
                    stem=base+"_"+arm
                    cost=matched[base] if arm=="pro-cost" else costs[stem]
                    raw=p/"raw"/((base+"_pro-cost-count") if arm=="pro-cost" else stem)
                    pool=raw.with_suffix(".pdbqt")
                    assert run.digest(pool)==cost["pool_sha256"],stem
                    k=dict(stem=stem,partition=("DYRK1A" if name in ("3ANR","7AKL") else "RBP4"),
                           query=query,receiver=rec,drug=name,group_molecule=query,seed=seed,arm=arm)
                    e=[]
                    for rank,block in enumerate(pool.read_text().split("ENDMDL")[:-1],1):
                        x=g.xyz(block,cell["serials"])
                        d=float(g.errors(x,cell["refs"]))
                        sc=float(next(l.split()[3] for l in block.splitlines() if l.startswith("REMARK VINA RESULT:")))
                        modes.append(dict(**k,rank=rank,rmsd_A=d,score=sc))
                        quality.append(dict(**k,rank=rank,**g.quality(x,src,mol,tree)))
                        e.append(d)
                    cells.append(dict(**k,n=len(e),top1_A=e[0] if e else "",best5_A=min(e[:5],default=""),
                                      oracle_A=min(e,default=""),regret_A=e[0]-min(e) if e else ""))
                    if arm not in run.ARMS:
                        continue
                    b=p/"dump"/stem
                    header,text=(b/"source.xyz").read_text().splitlines()
                    source_ids=list(map(int,header.split()[1:]))
                    x=np.array(text.split()[3:],float).reshape(-1,3)
                    order=np.argsort(source_ids)
                    x=x[order]
                    ids=np.array(source_ids)[order].tolist()
                    with np.load(b/"bank.npz") as bank:
                        poses=(x-x.mean(0))@bank["matched_rotation"].transpose(0,2,1)+bank["center"][:,None,:]
                        ds=g.errors(poses[:,[ids.index(s) for s in cell["serials"]]],cell["refs"])
                        ts=bank["triangle_id"]
                        res=bank["edge_residual_A2"]
                        fields={key:bank[key].item() for key in ("fitted_count","feasible_count","distance_evals","radius_A2")}
                    scores=run.old.csvrows(b/"scores.csv")
                    assert [int(v["index"]) for v in scores]==list(range(1024))
                    vals=np.array([float(v["score"]) for v in scores])
                    ranked=np.argsort(vals,kind="stable")
                    ranks=np.argsort(ranked)+1
                    idx=int(ranked[0])
                    inj=next(l.split() for l in (b/"inject.log").read_text().splitlines() if l.startswith("field_rank "))
                    assert int(inj[3])==idx and int(inj[11])==len(x)*1024
                    for j in range(len(ds)): cw.writerow((stem,j,ts[j],res[j],ds[j],vals[j],ranks[j]))
                    st,h,chain=g.inspect_fork(p,stem,cell,1,1)
                    stages.extend(st)
                    history.extend(h)
                    obs=(b/"stages.xyz").read_text().splitlines()
                    stage_ids=list(map(int,obs[0].split()[1:]))
                    first=next(line for line in obs if line.startswith("selected "))
                    actual=np.array(first.split()[3:],float).reshape(-1,3)
                    np.testing.assert_allclose(poses[idx,[ids.index(s) for s in stage_ids]],actual,atol=1e-9,rtol=0)
                    branches.append(dict(**k,index=idx,triangle_id=int(ts[idx]),n_le2=sum(ds<=2),generation_best_A=ds.min(),
                                first_le2_rank=min(ranks[ds<=2],default=""),selected_field_score=vals[idx],**fields,
                                second_A=next((v["rmsd_A"] for v in st if v["stage"]=="second"),""),
                                later_first_le2_step=next((v["index"] for v in st if v["stage"]=="offered" and v["scope"]==1 and v["index"]>1 and v["rmsd_A"]<=2),""),
                                **{k:v for k,v in chain.items() if k!="stem"}))
    for arm in (*run.ARMS,*native):
        for group in ("all","DYRK1A","RBP4",*run.CASES):
            rows=[r for r in cells if r["arm"]==arm and (group=="all" or group in (r["drug"],r["partition"]))]
            for unit in ("molecule","instance"):
                dist.append(historic.distribution(rows,arm,group,unit))
    keys=("top1_A","best5_A","oracle_A","regret_A")
    for row in cells:
        other=next(item for item in cells if item["drug"]==row["drug"] and item["seed"]==row["seed"] and item["arm"]=="edge")
        paired.append(dict(stem=row["stem"],reference=other["stem"],**{k:row[k]-other[k] if row["n"] and other["n"] else "" for k in keys}))
    for file,rows in (("cells",cells),("modes",modes),("chemistry",quality),("rank_stages",branches),
                    ("stages",stages),("residency",history),("distributions",dist),("paired",paired)):
        run.table(out/(file+".csv"),rows)
    run.old.write_json(out/"references.json",dict(raw_sha256=run.digest(run.BENCH/"cubic_grid/human_panel/depositions.tar.xz"),
                   serials={k:dict(serials=c["serials"],chemical_mappings=len(c["refs"])) for k,(c,*_) in refs.items()}))


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    for key in ("work","out"): cli.add_argument(key,type=lambda s:Path(s).resolve())
    a=cli.parse_args()
    evaluate(a.work,a.out)
