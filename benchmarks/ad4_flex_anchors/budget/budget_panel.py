"""Counterfactual reallocation of the cost of the AD4 anchor terms."""
import argparse
import csv
import itertools
import json
from pathlib import Path
import re
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from panel import GROUPS, SEEDS, ligand_only, measure, reference, run_timed, write_csv
from summarize import distribution

ARMS = ("old12000", "old13100", "old24000", "new12000")


def experiment(inputs, flex, experiment_folder, out, old, new):
    pools = out/"pools"
    pools.mkdir(parents=True, exist_ok=True)
    cases, modes = [], []
    for target, group in GROUPS.items():
        for nth,(q,r,seed) in enumerate(itertools.product(group,group,SEEDS)):
            meta, receiver = [json.loads((inputs / pdb / "preparation.json").read_text()) for pdb in (q,r)]
            expected=reference(meta,receiver)
            for arm in ARMS[nth%4:]+ARMS[:nth%4]:
                version,budget=arm[:3],int(arm[3:])
                stem=f"{q}-{r}-{seed}-{arm}"
                pose, log = pools/(stem+".pdbqt"), pools/(stem+".log")
                command=[new if version=="new" else old, "--maps",experiment_folder/"maps"/r/"receptor",
                         "--flex",flex/f"{r}_flex.pdbqt","--ligand",inputs/q/"ligand.pdbqt","--scoring","ad4",
                         "--seed",str(seed),"--exhaustiveness","1","--cpu","1","--num_modes","9",
                         "--energy_range","100","--max_evals",str(budget),"--out",pose,"--search_counters"]
                cost=run_timed(command,log)
                if cost["status"]:
                    raise RuntimeError(log.read_text())
                text=pose.read_text()
                if budget==12000 and pose.read_bytes()!=(experiment_folder/"pools"/f"{q}-{r}-{seed}-{version}.pdbqt").read_bytes():
                    raise ValueError(f"Original output differs for {stem}")
                sample=dict(target=target,query=q,receiver=r,seed=seed,arm=arm)
                errors=[]
                if "ENDMDL" in text:
                    only=pools/(stem+".ligand")
                    only.write_text(ligand_only(text))
                    errors=measure(only,meta["serial"],expected)
                    only.unlink()
                    scores=[float(line.split()[3]) for line in text.splitlines() if line.startswith("REMARK VINA RESULT:")]
                    assert len(scores)==len(errors)
                    for rank,(score,err) in enumerate(zip(scores,errors),1):
                        modes.append(dict(**sample,rank=rank,rmsd=err,score=score))
                trace=log.read_text()
                counters={key:int(re.search(r"\b"+key+r"=(\d+)",trace)[1]) for key in
                          ("hunt_evals","second_evals","explicit_evals")}
                cases.append(dict(**sample,budget=budget,**cost,**counters,reported=len(errors),
                                  top1=errors[0] if errors else "",best5=min(errors[:5]) if errors else "",
                                  best_reported=min(errors) if errors else "",regret=errors[0]-min(errors) if errors else ""))
    write_csv(out/"cases.csv",cases)
    write_csv(out/"modes.csv",modes)


def summarize(out):
    cases=list(csv.DictReader((out/"cases.csv").open()))
    rows,pairs,intervals,costs,ecdf=[],[],[],[],[]
    for target in ("all",*GROUPS):
        block=[x for x in cases if target=="all" or x["target"]==target]
        for arm in ARMS:
            sub=[x for x in block if x["arm"]==arm]
            costs.append(dict(target=target,arm=arm,n=len(sub),success=sum(bool(s["top1"]) for s in sub),
                              wall_s=sum(float(s["wall_s"]) for s in sub),cpu_s=sum(float(s["cpu_s"]) for s in sub),
                              median_rss_mib=np.median([float(s["peak_kib"])/1024 for s in sub]),
                              total_evals=sum(int(s[k]) for s in sub for k in ("hunt_evals","second_evals","explicit_evals")),
                              pose_count=sum(int(s["reported"]) for s in sub)))
        for metric in ("top1","best5","best_reported","regret"):
            arms={arm:{(s["query"],s["receiver"],int(s["seed"])):float(s[metric]) for s in block
                       if s["arm"]==arm and s[metric]!=""} for arm in ARMS}
            for arm,values in arms.items():
                observed=list(values.values())
                rows.append(dict(target=target,metric=metric,arm=arm,**distribution(observed)))
                support,counts=np.unique(observed,return_counts=True)
                for x,n in zip(support,np.cumsum(counts)):
                    ecdf.append(dict(target=target,metric=metric,arm=arm,rmsd=x,count_le=n,total=len(observed),
                                     ecdf=n/len(observed)))
            for test,base in (("new12000","old13100"),("new12000","old24000"),("new12000","old12000"),
                              ("old13100","old12000"),("old24000","old12000")):
                label=test+"-minus-"+base
                diff={k:arms[test][k]-arms[base][k] for k in sorted(arms[test].keys() & arms[base].keys())}
                rows.append(dict(target=target,metric=metric,arm=label,**distribution(list(diff.values()))))
                for (q,r,seed),value in diff.items():
                    pairs.append(dict(target=target,metric=metric,comparison=label,query=q,receiver=r,seed=seed,delta=value))
                sm={s:np.mean([v for k,v in diff.items() if k[2]==s]) for s in SEEDS}
                boot=[np.mean([sm[s] for s in draws]) for draws in itertools.product(SEEDS,repeat=len(SEEDS))]
                intervals.append(dict(target=target,metric=metric,comparison=label,paired_delta=np.mean(list(diff.values())),
                                      lower=np.percentile(boot,2.5),upper=np.percentile(boot,97.5)))
    write_csv(out/"distributions.csv",rows)
    write_csv(out/"paired.csv",pairs)
    write_csv(out/"intervals.csv",intervals)
    write_csv(out/"costs.csv",costs)
    write_csv(out/"ecdf.csv",ecdf)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("inputs","flex","experiment","out","old","new"):
        p.add_argument(arg,type=Path)
    a=p.parse_args()
    experiment(*(getattr(a,k).resolve() for k in ("inputs","flex","experiment","out","old","new")))
    summarize(a.out)
