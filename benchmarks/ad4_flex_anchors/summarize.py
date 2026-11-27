"""Derive paired distributions and fixed-pool reranking from the archived raw rows."""
import argparse
from collections import defaultdict
import csv
import itertools
from pathlib import Path

import numpy as np
from panel import SEEDS, write_csv


def distribution(x):
    return dict(zip(("n","mean","min","q25","median","q75","p90","p95","max"),
                    [len(x), np.mean(x), *np.percentile(x,[0,25,50,75,90,95,100])]))


def load(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder",type=Path)
    a=p.parse_args().folder
    cases=load(a/"cases.csv")
    rows=[]
    deltas=[]
    uncertainty=[]
    for target in ("all","FABP4","RBP4"):
        for metric in ("top1","best5","best_reported","regret"):
            by_arm={}
            for arm in ("old","new","official"):
                by_arm[arm]={(r["query"],r["receiver"],int(r["seed"])):float(r[metric]) for r in cases
                             if (target=="all" or r["target"]==target) and r["arm"]==arm and r[metric]!=""}
                rows.append(dict(target=target,metric=metric,arm=arm,**distribution(list(by_arm[arm].values()))))
            common=sorted(by_arm["new"].keys() & by_arm["old"].keys())
            diff={k:by_arm["new"][k]-by_arm["old"][k] for k in common}
            for k,d in diff.items():
                deltas.append(dict(target=target,metric=metric,query=k[0],receiver=k[1],seed=k[2],delta=d))
            rows.append(dict(target=target,metric=metric,arm="new-minus-old",**distribution(list(diff.values()))))
            # Resample entire seed blocks simultaneously across the fixed receptor/query grid.
            seedmeans={s:np.mean([d for k,d in diff.items() if k[2]==s]) for s in SEEDS}
            boot=[np.mean([seedmeans[s] for s in draws]) for draws in itertools.product(SEEDS,repeat=len(SEEDS))]
            uncertainty.append(dict(target=target,metric=metric,independent_seed_blocks=len(SEEDS),
                                    paired_mean=np.mean(list(diff.values())),conditional_lo=np.percentile(boot,2.5),
                                    conditional_hi=np.percentile(boot,97.5)))
    write_csv(a/"distributions.csv",rows)
    write_csv(a/"paired.csv",deltas)
    write_csv(a/"conditional_intervals.csv",uncertainty)
    pools=defaultdict(list)
    for r in load(a/"fixed.csv"):
        pools[r["query"],r["receiver"],r["seed"]].append(r)
    fixed=[]
    for (q,r,seed), pool in pools.items():
        sel = [min(pool,key=lambda row:(float(row[arm+"_score"]),int(row["rank"]))) for arm in ("old","new")]
        fixed.append(dict(query=q,receiver=r,seed=seed,n=len(pool),old_rank=sel[0]["rank"],new_rank=sel[1]["rank"],
                          old_rmsd=sel[0]["rmsd"],new_rmsd=sel[1]["rmsd"],
                          delta_rmsd=float(sel[1]["rmsd"])-float(sel[0]["rmsd"])))
    write_csv(a/"fixed_choices.csv",fixed)

if __name__ == "__main__":
    main()
