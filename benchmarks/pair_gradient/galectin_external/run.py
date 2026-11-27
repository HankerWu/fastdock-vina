"""Score the frozen external Gal3 series by native top pose at one seed."""
import argparse
import csv
import hashlib
import importlib.util
import itertools
import json
import os
import random
from pathlib import Path
from statistics import median

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("gal_time", HERE.parent / "galectin_time/run.py")
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("--parent",type=Path,required=True)
    p.add_argument("--linear",type=Path,required=True)
    p.add_argument("--cpu",type=int,default=8)
    a=p.parse_args()
    os.sched_setaffinity(0,{a.cpu})
    w=a.work.resolve()
    logs=w/"docking"
    (logs/"raw").mkdir(parents=True,exist_ok=True)
    with (w/"manifest.csv").open() as f:
        rows=list(csv.DictReader(f))
    bins={arm:getattr(a,arm).resolve() for arm in ("parent","linear")}
    with (HERE.parent/"galectin_time/budgets.csv").open() as f:
        cal=list(csv.DictReader(f))
    caps={rec:int(median(int(r["calibrated_cap"]) for r in cal if r["receptor"]==rec))
          for rec in ("1KJL","3T1L")}
    cells=list(itertools.product(caps,[int(r["compound"]) for r in rows]))
    random.Random(214039).shuffle(cells)
    primary=[]
    for i,(rec,comp) in enumerate(cells):
        for arm in (("parent","linear") if i%2 else ("linear","parent")):
            primary.append(prior.dock(logs,w/"prepared",bins,rec,comp,214039,arm,12000))
    prior.save(w/"fixed.csv",primary)
    def total(arm,keys):
        return sum(sum(r[k] for k in keys) for r in primary if r["arm"]==arm)
    ratios=[total("linear",k) / total("parent",k) for k in
            (["elapsed_ns"],["hunt_evals","second_evals","explicit_evals"])]
    cache={(r["receptor"],r["compound"],r["arm"],r["cap"]):r for r in primary}
    def frozen(rec,comp,cap):
        key=rec,comp,"linear",cap
        if key not in cache:
            cache[key]=prior.dock(logs,w/"prepared",bins,rec,comp,214039,"linear",cap)
        return cache[key]
    transported=[frozen(rec,comp,caps[rec]) for rec,comp in cells]
    prior.save(w/"transported.csv",transported)
    if any(abs(r-1)>.05 for r in ratios):
        cal_rows=[]; cap_rows=[]
        for i,(rec,comp) in enumerate(cells):
            arms=[("parent",12000),("linear",3000),("linear",12000)]
            times={}
            pos=i%3
            for arm,cap in arms[pos:]+arms[:pos]:
                r=prior.dock(logs,w/"prepared",bins,rec,comp,214033,arm,cap)
                times[arm,cap]=r["elapsed_ns"]
                cal_rows.append(r)
            low,high,target=times["linear",3000],times["linear",12000],times["parent",12000]
            cap=12000 if high<=target else 3000 if low>=target or high<=low else 3000+9000*(target-low)//(high-low)
            cap_rows.append({"receptor":rec,"compound":comp,"cap":cap})
        prior.save(w/"calibration.csv",cal_rows); prior.save(w/"budgets.csv",cap_rows)
        prior.save(w/"matched.csv",[frozen(r["receptor"],r["compound"],r["cap"]) for r in cap_rows])
    (w/"run.json").write_text(json.dumps({"seed":214039,"transported_caps":caps,"cpu":a.cpu,
        "fixed_linear_parent_time_joint_ratios":ratios,"series_specific_calibration":any(abs(r-1)>.05 for r in ratios),
        "binary_sha256":{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in bins.items()},
        "manifest_sha256":hashlib.sha256((w/"manifest.csv").read_bytes()).hexdigest()},indent=2)+"\n")


if __name__=="__main__":
    main()
