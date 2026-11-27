"""Conditional extra native trajectories at the saturated cap, per SENSITIVITY.md."""
import argparse
import os
import math
from pathlib import Path
import resource
import subprocess
import time

import run_panel as panel


def native(p,c,base,seed,n,binary,count=False):
    stem=f"{base}_pro-e{n}" + ("-count" if count else "")
    pool=p/"raw"/(stem+".pdbqt")
    args=panel.parent.native(c,seed,pool)
    args[args.index("--exhaustiveness")+1]=n
    args[args.index("--max_evals")+1]=800000
    log=pool.with_suffix(".log")
    start,prev=time.perf_counter(),resource.getrusage(resource.RUSAGE_CHILDREN)
    env={k:v for k,v in os.environ.items() if not k.startswith("FD_")}
    with log.open("w") as f: subprocess.run(list(map(str,[binary,*args])),check=True,stdout=f,stderr=f,env=env)
    end,cpu=time.perf_counter(),resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(stem=stem,wall_s=end-start,cpu_s=cpu.ru_utime+cpu.ru_stime-prev.ru_utime-prev.ru_stime,
                pool_sha256=panel.digest(pool),**panel.parent.counts(log),
                **(panel.old.interpolation(log) if count else {}))


def predict(p,release,observer):
    prior={r["case"]:r for r in panel.old.csvrows(p/"cost_selected.csv")}
    observed,checks,chosen=[],[],[]
    for name,c in panel.old.configs(p).items():
        for seed in panel.CASES[name][2]:
            base=f"{name}_{seed}"
            s=prior[base]
            t=float(s["target_cpu_s"])
            orig=int(s["cap"])
            options={2:dict(stem=base+"_pro-cost-count",**{k:s[k] for k in ("wall_s","cpu_s","pool_sha256")})}
            n=2
            if not int(s["bracketed"]):
                assert orig==800000 and int(s["cap_unused"])
                while float(options[n]["cpu_s"])<t:
                    n+=1
                    row=native(p,c,base,seed,n,release)
                    options[n]=row
                    observed.append(row)
                    panel.table(p/"extra_trajectory_runs.csv",observed)
            best=min(options,key=lambda k:(abs(math.log(float(options[k]["cpu_s"])/t)),k))
            win=options[best]
            if best>2:
                check=native(p,c,base,seed,best,observer,count=True)
                assert check["pool_sha256"]==win["pool_sha256"]
                checks.append(check)
                panel.table(p/"extra_trajectory_count.csv",checks)
            chosen.append(dict(case=base,selected_stem=win["stem"],trajectories=best,cap=orig,
                               target_cpu_s=t,cpu_s=float(win["cpu_s"]),wall_s=float(win["wall_s"]),
                               cpu_ratio=float(win["cpu_s"])/t,pool_sha256=win["pool_sha256"]))
            panel.table(p/"extra_trajectory_selected.csv",chosen)


def evaluate(p,out):
    # Load labels only after every duration-based choice is recorded.
    from evaluate import g
    refs=g.reference_inputs(p)
    values,rows,total=[],[],[]
    counters={r["stem"]:r for fname in ("ladder_counters","extra_trajectory_count")
              for r in panel.old.csvrows(p/(fname+".csv"))}
    for r in panel.old.csvrows(p/"extra_trajectory_selected.csv"):
        name,_=r["case"].split("_")
        cell,mol,src,tree=refs[name]
        pool=p/"raw"/(r["selected_stem"]+".pdbqt")
        assert panel.digest(pool)==r["pool_sha256"]
        counter_stem=r["selected_stem"] if r["selected_stem"].endswith("-count") else r["selected_stem"]+"-count"
        total.append(dict(case=r["case"],wall_s=r["wall_s"],cpu_s=r["cpu_s"],**{k:int(counters[counter_stem][k])
                          for k in (*panel.old.COUNT_KEYS,"gradient_fields","scalar_fields")}))
        errs=[]
        for rank,block in enumerate(pool.read_text().split("ENDMDL")[:-1],1):
            x=g.xyz(block,cell["serials"])
            error=float(g.errors(x,cell["refs"]))
            errs.append(error)
            values.append(dict(case=r["case"],rank=rank,rmsd_A=error,**g.quality(x,src,mol,tree)))
        rows.append(dict(**r,n=len(errs),top1_A=errs[0] if errs else "",best5_A=min(errs[:5],default=""),
                         oracle_A=min(errs,default=""),regret_A=errs[0]-min(errs) if errs else ""))
    panel.table(out/"extra_trajectory_cost.csv",total)
    panel.table(out/"extra_trajectory_cells.csv",rows)
    panel.table(out/"extra_trajectory_modes.csv",values)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("action",choices=("predict","evaluate"))
    cli.add_argument("work",type=lambda s:Path(s).resolve())
    cli.add_argument("--core",type=int)
    for key in ("production","pro_counter","out"): cli.add_argument("--"+key,type=lambda s:Path(s).resolve())
    a=cli.parse_args()
    if a.action=="predict":
        os.sched_setaffinity(0,{a.core})
        predict(a.work,a.production,a.pro_counter)
    else: evaluate(a.work,a.out)
