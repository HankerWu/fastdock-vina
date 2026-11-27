"""Validate prediction replays, paired source identities, and cost counts from native logs."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import sys

import panel


def run(p, out):
    traced = {r["stem"]:r for r in panel.old.csvrows(p/"run_cost.csv")}
    timed = {r["stem"]:r for r in panel.old.csvrows(p/"live_cost.csv")}
    ctl = panel.old.csvrows(p/"control_cost.csv")
    release = {r["stem"].removesuffix("-production"):r for r in panel.old.csvrows(p/"production_cost.csv")}
    cell_stats, sources = [], []
    def row(stem, arm, work, measured, site=(0,0,0,0)):
        assert work["pool_sha256"] == measured["pool_sha256"]
        shape,screen,nearest,burial = site
        return dict(stem=stem,arm=arm,wall_s=float(measured["wall_s"]),cpu_s=float(measured["cpu_s"]),
                    site_nearest=nearest,site_burial=burial,shape_queries=shape,rank_field_queries=screen,
                    **{k:int(work[k]) for k in ("hunt_calls","hunt_evals","second_calls","second_evals","explicit_calls",
                                                "explicit_evals","late_prescore","scalar_fields","gradient_fields")})
    for name,cfg in panel.old.configs(p).items():
        for seed in panel.CASES[name][2]:
            base=f"{name}_{seed}"
            for suffix in ("source.xyz","source.q"):
                hashes = [panel.digest(p/"dump"/(base+"_"+arm)/suffix) for arm in panel.ARMS]
                assert len(set(hashes))==1,(base,suffix)
                sources.append(dict(case=base,file=suffix,sha256=hashes[0]))
            for arm in panel.ARMS:
                stem=base+"_"+arm
                cli = (p/"timed"/(stem+".log")).read_text()
                capture = (p/"raw"/(stem+".log")).read_text()
                def diagnostic(text):
                    return [l for l in text.splitlines() if l.startswith(("Search trajectory:","Search counters:",
                            "Hunt chain:","site_nearest ","field_rank ","FIELD_COUNT "))]
                assert diagnostic(cli) == diagnostic(capture), stem
                site = next(l.split() for l in cli.splitlines() if l.startswith("site_nearest "))
                screen = int(next(l.split() for l in cli.splitlines() if l.startswith("field_rank "))[11])
                cell_stats.append(row(stem,arm,traced[stem],timed[stem],(int(site[5]),screen,int(site[1]),int(site[3]))))
    for r in ctl:
        stem=r["stem"]
        arm=stem.rsplit("_",1)[1]
        measured = r
        if arm.startswith("pro"):
            measured = timed[stem] if arm=="pro" else release[stem]
        cell_stats.append(row(stem,arm,r,measured))
    panel.table(out/"source_identity.csv",sources)
    panel.table(out/"cost_cells.csv",cell_stats)
    total=[]
    for arm in sorted({r["arm"] for r in cell_stats}):
        group=[r for r in cell_stats if r["arm"]==arm]
        total.append(dict(arm=arm,cells=len(group),**{k:sum(r[k] for r in group) for k in group[0] if k not in ("stem","arm")}))
    panel.table(out/"cost_total.csv",total)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ("work","out"): parser.add_argument(k,type=Path)
    args=parser.parse_args()
    run(args.work.resolve(),args.out.resolve())
