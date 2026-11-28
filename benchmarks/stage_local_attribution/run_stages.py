"""Execute the pre-registered stage cross and independent replay of every stable-path source."""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "bfgs_source_replay"))
from run_replay import run_panel, digest
from analyze import native

CASES = {("7AKL","3ANR",85237), ("5TBM","4GHI",85231), ("1FKH","1FKG",85237)}
MODES = ("00", "10", "01", "11")


def cases(work):
    for job in run_panel.panel(work / "fixtures"):
        cohort, inputs, q, r, ref, seed, cap = job
        if cohort == "reciprocal" and (q["pdb"], r["pdb"], seed) in CASES:
            yield f"{q['pdb']}-{r['pdb']}-{seed}", job


def replay_inputs(path):
    trajectories, _ = native(path)
    for trajectory in trajectories:
        for events in trajectory:
            yield "H", events["P"][0][3], events["H"][0]
            if events["D"]:
                yield "D", events["H"][0][3], events["D"][0]


def main(work, binary):
    work.mkdir(parents=True, exist_ok=False)
    archive = list(csv.DictReader((HERE.parent/"search_trial_panel/results/cells.csv").open()))
    costs = []
    def record(row):
        costs.append(row)
        (work/"costs.json").write_text(json.dumps(costs,indent=2)+"\n")
    for case, job in cases(work):
        cohort, inputs, q, r, ref, seed, cap = job
        casepath = work / case
        casepath.mkdir()
        for mode in MODES:
            arm = "stable" if mode != "11" else "candidate"
            dest = casepath / mode
            dest.mkdir()
            trace = dest / "free.native"
            os.environ["FD_BFGS_TRACE"] = str(trace)
            os.environ["FD_STAGE_MODE"] = mode
            os.environ.pop("FD_BFGS_SHADOW",None)
            stats, ranks = run_panel.dock((arm, binary, dest, job))
            record(dict(case=case, mode=mode, phase="free", binary_sha256=digest(binary), **stats))
            assert stats["status"] == "ok"
            if mode in ("00","11"):
                row = next(row for row in archive if row["key"] == stats["key"])
                assert row["pool_sha256"] == stats["pool_sha256"], (row, stats)
            run_panel.csvfile(dest/"published.csv",ranks)
        inputs_path = casepath/"stable.states"
        baseline = list(replay_inputs(casepath/"00/free.native"))
        inputs_path.write_text("".join(f"{event} {source}\n" for event,source,_ in baseline))
        for mode in ("00","11"):
            trace = casepath / mode / "locked.trace"
            log = trace.with_suffix(".log")
            os.environ.update(FD_STAGE_MODE=mode,FD_BFGS_SHADOW=str(inputs_path),FD_BFGS_TRACE=str(trace))
            cmd = [str(binary), "--ligand", str(inputs/q["pdb"]/"ligand.pdbqt"), "--receptor",
                   str(inputs/r["pdb"]/"receptor.pdbqt"), "--cpu","1", "--seed",str(seed),
                   "--max_evals",str(cap),"--exhaustiveness","3","--num_modes","9","--energy_range","100",
                   "--spacing","0.375","--out",str(trace.with_suffix(".unused"))]
            for axis,center in zip("xyz",r["center"]):
                cmd += ["--size_"+axis,"24","--center_"+axis,str(center)]
            stats = run_panel.execute(cmd,log,180)
            record(dict(case=case,mode=mode,phase="locked",binary_sha256=digest(binary),cmd=cmd,**stats))
            assert stats["status"] == "ok"
            rows = [line.split("\t") for line in trace.read_text().splitlines()]
            assert len(rows) == len(baseline)
            for i, ((event,source,expected),rec) in enumerate(zip(baseline, rows)):
                assert rec[0] == event and int(rec[1]) == i
                assert mode != "00" or expected[2:-1] == rec[2:-1], (case, i)
            hunt = [rec for rec in rows if rec[0]=="H"]
            print(case,mode,len(hunt),sum(int(rec[4]) for rec in hunt),
                  sum(int(rec[13]) for rec in hunt), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("binary",type=Path)
    a=p.parse_args()
    main(a.work.resolve(),a.binary.resolve())
