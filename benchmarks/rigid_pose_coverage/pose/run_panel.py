"""Prospective predictions. This module does not load query crystals or evaluation labels."""
import argparse
import math
import os
import re
from pathlib import Path
import resource
import subprocess
import sys
import time
import importlib.util

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
spec = importlib.util.spec_from_file_location("fragment_panel", BENCH/"fragment_anchor/pose/panel.py")
parent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parent)
old = parent.old
table, digest, invoke = parent.table, parent.digest, parent.invoke

CASES = {
    "3ANR": ("3ANR", "7AKL", (96107, 96131)),
    "7AKL": ("7AKL", "3ANR", (96107, 96131)),
    "3FMZ-A": ("3FMZ", "1RBP", (96107,)),
    "3FMZ-HA": ("3FMZ", "1RBP", (96107,)),
}
parent.CASES = old.CASES = CASES
ARMS = ("edge", "feasible_edge", "coverage", "uniform", "uniform_work")


def live(p, c, seed, arm, binary, trace):
    stem = f"{c['ligand'].stem.removesuffix('-ligand')}_{seed}_{arm}"
    raw, dump = p/"raw"/stem, p/"dump"/stem
    out = (raw if trace else p/"timed"/stem).with_suffix(".pdbqt")
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, HERE.parent/"engine/live.py", "--engine", binary, "--arm", arm]
    if trace:
        cmd += ["--dump", dump]
    cmd += ["--", *parent.native(c, seed, out)]
    log = out.with_suffix(".log")
    env = {k:v for k,v in os.environ.items() if not k.startswith("FD_")}
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    with log.open("w") as stream:
        tic, cpu0 = time.perf_counter(), resource.getrusage(resource.RUSAGE_CHILDREN)
        subprocess.run(list(map(str,cmd)), check=True, env=env, stdout=stream, stderr=subprocess.STDOUT)
        toc, cpu1 = time.perf_counter(), resource.getrusage(resource.RUSAGE_CHILDREN)
    if trace:
        # The shared native-lineage reader opens flat stem files.
        for ext, file in (("stages","stages.xyz"),("history","history"),("inject","inject.log")):
            raw.with_suffix("."+ext).symlink_to(Path("..")/"dump"/stem/file)
    else:
        assert digest(out) == digest(raw.with_suffix(".pdbqt")), stem
    return dict(stem=stem, wall_s=toc-tic,
                cpu_s=cpu1.ru_utime+cpu1.ru_stime-cpu0.ru_utime-cpu0.ru_stime,
                pool_sha256=digest(out), **parent.counts(log), **old.interpolation(log))


def dock(p, engine, pro, trace):
    rows=[]
    for name,c in old.configs(p).items():
        for seed in CASES[name][2]:
            arms=ARMS if trace else (*ARMS,"pro")
            shift=seed%len(arms)
            for arm in arms[shift:]+arms[:shift]:
                if arm=="pro":
                    cache=[]
                    stem=f"{name}_{seed}_pro"
                    out=invoke(p,cache,stem+"-paired-clock",pro,c,seed)
                    assert digest(out)==digest(p/"raw"/(stem+".pdbqt"))
                    item=cache[0]|dict(stem=stem,scalar_fields="",gradient_fields="")
                else:
                    item=live(p,c,seed,arm,engine,trace)
                rows.append(item)
                table(p/("run_cost.csv" if trace else "live_cost.csv"), rows)
                print(item["stem"], flush=True)


def ladder(p, production, observer):
    """Blindly bracket coverage CPU using native searches and verify counters on the chosen run."""
    measured={r["stem"]:r for r in old.csvrows(p/"live_cost.csv")}
    previous={r["stem"].removesuffix("-production"):r for r in old.csvrows(p/"production_cost.csv")}
    summaries, runs, counters=[], [], []
    for name,c in old.configs(p).items():
        for seed in CASES[name][2]:
            base=f"{name}_{seed}"
            target=float(measured[base+"_coverage"]["cpu_s"])
            options={25000: measured[base+"_pro"],50000: previous[base+"_pro-high"]}
            cap=100000
            unused=False
            while cap<=400000 or (float(options[max(options)]["cpu_s"])<target and not unused):
                invoke(p,runs,f"{base}_pro-{cap}",production,c,seed,cap)
                options[cap]=runs[-1]
                log=(p/"raw"/(runs[-1]["stem"]+".log")).read_text()
                work=re.findall(r"^Search trajectory: .*?hunt_evals=(\d+) .*?second_evals=(\d+)",log,re.M)
                assert len(work)==2
                # Both trajectories finished below the cap; increasing it cannot add any evaluations.
                unused=all(int(h)+int(s)<cap for h,s in work)
                cap*=2
            select=min(options,key=lambda k:(abs(math.log(float(options[k]["cpu_s"])/target)),k))
            chosen=options[select]
            pool=invoke(p,counters,base+"_pro-cost-count",observer,c,seed,select)
            assert digest(pool)==chosen["pool_sha256"]
            counters[-1].update(old.interpolation(pool.with_suffix(".log")))
            summaries.append(dict(case=base,arm="pro-cost",cap=select,wall_s=chosen["wall_s"],cpu_s=chosen["cpu_s"],
                                  target_cpu_s=target,bracketed=int(min(float(v["cpu_s"]) for v in options.values())<=target<=max(float(v["cpu_s"]) for v in options.values())),
                                  max_cap=max(options),cap_unused=int(unused),pool_sha256=chosen["pool_sha256"]))
            table(p/"ladder.csv",runs)
            table(p/"ladder_counters.csv",counters)
            table(p/"cost_selected.csv",summaries)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("action",choices=("prepare","controls","dock","timing","ladder"))
    cli.add_argument("work",type=lambda s:Path(s).resolve())
    cli.add_argument("--core",type=int,required=True)
    for key in ("observer","counter","official","production","pro_counter"):
        cli.add_argument("--"+key,type=lambda s:Path(s).resolve())
    a=cli.parse_args()
    os.sched_setaffinity(0,{a.core})
    if a.action=="prepare": old.prepare(a.work)
    elif a.action=="controls": parent.controls(a.work,a.production,a.pro_counter,a.official,a.counter)
    elif a.action=="ladder": ladder(a.work,a.production,a.pro_counter)
    else: dock(a.work,a.observer,a.production,a.action=="dock")
