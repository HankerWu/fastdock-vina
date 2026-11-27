"""Blinded pose inputs and runs. Reference structures are read only by evaluate.py."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
RANK = BENCH / "cavity_field_rank"
spec = importlib.util.spec_from_file_location("old_field_panel", RANK / "pose/run.py")
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
from pose_runs import counts, digest, invoke, table

CASES = {
    **old.CASES,
    "3FMZ-HA": ("3FMZ", "1RBP", (91813, 93241, 94709)),
    "3FMZ-A": ("3FMZ", "1RBP", (91813, 93241, 94709)),
    "2NNQ": ("2NNQ", "3FR2", (94709, 94727)),
    "3FR2": ("3FR2", "2NNQ", (94709, 94727)),
}
old.CASES = CASES
ARMS = ("whole", "fragment", "sham")


def partition(name, seed):
    if name in ("2NNQ", "3FR2"):
        return "transfer_FABP4"
    if seed == 94709:
        return "transfer_A1120"
    return "previous10"


def controls(p, pro, observer, official, counter):
    old.controls(p, observer, official, counter)
    timed = []
    for name, c in old.configs(p).items():
        for seed in CASES[name][2]:
            for arm, cap in (("pro", 25000), ("pro-high", 50000)):
                stem = f"{name}_{seed}_{arm}"
                out = invoke(p, timed, stem+"-production", pro, c, seed, cap)
                assert out.read_bytes() == (p/"raw"/(stem+".pdbqt")).read_bytes()
    table(p/"production_cost.csv",timed)
    (p/"production.sha256").write_text(digest(pro)+"\n")


def native(c, seed, out):
    args = ["--ligand", c["ligand"], "--ligand_bonds", c["bonds"], "--receptor", c["receptor"],
            "--cpu", 1, "--seed", seed, "--max_evals", 25000, "--exhaustiveness", 2,
            "--num_modes", 9, "--energy_range", 100, "--pose_selection", "late", "--search_counters", "--out", out]
    for axis, center in zip("xyz",c["center"],strict=True):
        args.extend(["--center_"+axis, center, "--size_"+axis, c["size"]])
    return args


def dock(p, engine, trace=True, production=None):
    costs = []
    expected = {r["stem"]:r["pool_sha256"] for r in old.csvrows(RANK/"pose/results/dock_cost.csv")}
    for name,c in old.configs(p).items():
        for j,seed in enumerate(CASES[name][2]):
            arms = (*ARMS,"pro") if production is not None else ARMS
            for arm in arms[j%len(arms):]+arms[:j%len(arms)]:
                if arm == "pro":
                    timed = []
                    tag = f"{name}_{seed}_pro"
                    out = invoke(p,timed,tag+"-paired-clock",production,c,seed)
                    assert out.read_bytes() == (p/"raw"/(tag+".pdbqt")).read_bytes()
                    costs.append(dict(stem=tag,arm=arm,seed=seed,**{k:v for k,v in timed[0].items() if k!="stem"},
                                      scalar_fields="",gradient_fields=""))
                    table(p/"live_cost.csv",costs)
                    continue
                stem = f"{name}_{seed}_{arm}"
                raw = p/"raw"/stem
                root = p/"dump"/stem
                parent = p/("raw" if trace else "timed")
                parent.mkdir(exist_ok=True)
                out = parent/(stem+".pdbqt")
                cmd = [sys.executable, HERE.parent/"engine/live.py", "--engine",engine, "--arm", arm]
                if trace:
                    root.parent.mkdir(exist_ok=True)
                    cmd.extend(["--dump", root])
                cmd.extend(["--", *native(c,seed,out)])
                log = out.with_suffix(".log")
                env = {k:v for k,v in os.environ.items() if not k.startswith("FD_")}
                env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
                with log.open("w") as stream:
                    t0 = time.perf_counter()
                    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
                    subprocess.run(list(map(str,cmd)), stdout=stream, stderr=subprocess.STDOUT, env=env, check=True)
                    t1 = time.perf_counter()
                    end = resource.getrusage(resource.RUSAGE_CHILDREN)
                if trace:
                    # Shared lineage reader accepts the historical flat names.
                    for ext, orig in (("stages","stages.xyz"),("inject","inject.log"),("history","history")):
                        dest = raw.with_suffix("."+ext)
                        dest.unlink(missing_ok=True)
                        dest.symlink_to(Path("..")/"dump"/stem/orig)
                else:
                    assert digest(out) == digest(raw.with_suffix(".pdbqt")), stem
                if arm == "whole" and partition(name,seed) == "previous10":
                    assert digest(out) == expected[f"{name}_{seed}_true"], stem
                counts1 = old.interpolation(log)
                costs.append(dict(stem=stem,arm=arm,seed=seed,wall_s=t1-t0,
                                  cpu_s=end.ru_utime+end.ru_stime-usage.ru_utime-usage.ru_stime,
                                  pool_sha256=digest(out),**counts(log),**counts1))
                table(p/("run_cost.csv" if trace else "live_cost.csv"),costs)
                print(stem,flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare","controls","dock","timing"))
    parser.add_argument("work", type=lambda s:Path(s).resolve())
    for key in ("pro","observer","official","counter","core"):
        parser.add_argument("--"+key,type=int if key=="core" else lambda s:Path(s).resolve())
    args = parser.parse_args()
    if args.core is not None:
        os.sched_setaffinity(0,{args.core})
    if args.action == "prepare":
        old.prepare(args.work)
    if args.action == "controls":
        controls(args.work,args.pro,args.observer,args.official,args.counter)
    if args.action in ("dock","timing"):
        dock(args.work,args.observer,trace=args.action=="dock",production=args.pro if args.action=="timing" else None)
