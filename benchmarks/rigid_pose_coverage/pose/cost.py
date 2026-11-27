"""Read timings, counters and hashes; require identical sources and released pools on replay."""
import argparse
from pathlib import Path

import numpy as np
import run_panel as panel


def rows_by(p,file):
    return {r["stem"]:r for r in panel.old.csvrows(p/(file+".csv"))}


def audit(p,out):
    traced=rows_by(p,"run_cost")
    timed=rows_by(p,"live_cost")
    baseline=rows_by(p,"control_cost")
    release=rows_by(p,"production_cost")
    selected={r["case"]:r for r in panel.old.csvrows(p/"cost_selected.csv")}
    counts=rows_by(p,"ladder_counters")
    rows,sources=[],[]
    def add(stem,arm,measured,queries,geometry=None):
        assert queries["pool_sha256"]==measured["pool_sha256"]
        g = dict(site_nearest=0,site_burial=0,shape_queries=0,rank_field_queries=0,
                 fits=0,distance_evals=0)
        if geometry: g.update(geometry)
        rows.append(dict(stem=stem,arm=arm,wall_s=float(measured["wall_s"]),cpu_s=float(measured["cpu_s"]),
                         **g,**{key:int(queries[key]) for key in (*panel.old.COUNT_KEYS,"scalar_fields","gradient_fields")}))
    for name,(_,_,seeds) in panel.CASES.items():
        for seed in seeds:
            base=f"{name}_{seed}"
            for suffix in ("source.xyz","source.q"):
                hashes={panel.digest(p/"dump"/(base+"_"+arm)/suffix) for arm in panel.ARMS}
                assert len(hashes)==1,(base,suffix)
                sources.append(dict(case=base,filename=suffix,sha256=hashes.pop()))
            for file in ("history","inject.log","scores.csv","stages.xyz"):
                assert (p/"dump"/(base+"_uniform")/file).read_bytes()==(p/"dump"/(base+"_uniform_work")/file).read_bytes()
            for arm in panel.ARMS:
                stem=base+"_"+arm
                log=(p/"timed"/(stem+".log")).read_text()
                tags=("Hunt chain:","Search counters:","Search trajectory:","site_nearest ","field_rank ","FIELD_COUNT ")
                diagnostic=lambda s:[line for line in s.splitlines() if line.startswith(tags)]
                assert diagnostic(log)==diagnostic((p/"raw"/(stem+".log")).read_text())
                site=next(l.split() for l in log.splitlines() if l.startswith("site_nearest "))
                field=next(l.split() for l in log.splitlines() if l.startswith("field_rank "))
                with np.load(p/"dump"/stem/"bank.npz") as b:
                    add(stem,arm,timed[stem],traced[stem],dict(site_nearest=int(site[1]),site_burial=int(site[3]),
                        shape_queries=int(site[5]),rank_field_queries=int(field[11]),fits=int(b["fitted_count"]),
                        distance_evals=int(b["distance_evals"])))
            for arm in ("pro","pro-high","official","official-high","pro-cost"):
                stem=base+"_"+arm
                if arm=="pro-cost":
                    add(stem,arm,selected[base],counts[base+"_pro-cost-count"])
                else:
                    measured=timed[stem] if arm=="pro" else (release[stem+"-production"] if arm=="pro-high" else baseline[stem])
                    add(stem,arm,measured,baseline[stem])
    total=[]
    for arm in (*panel.ARMS,"pro","pro-high","official","official-high","pro-cost"):
        g=[r for r in rows if r["arm"]==arm]
        total.append(dict(arm=arm,n=len(g),**{k:sum(r[k] for r in g) for k in g[0] if k not in ("arm","stem")}))
    for f,data in (("source_identity",sources),("cost_cells",rows),("cost_total",total)):
        panel.table(out/(f+".csv"),data)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    for key in ("work","out"): cli.add_argument(key,type=lambda s:Path(s).resolve())
    a=cli.parse_args()
    audit(a.work,a.out)
