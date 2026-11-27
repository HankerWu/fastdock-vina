"""Compare one-pass docking against all frozen-bank published coordinates and counters."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
SHAPE_ADAPTER = HERE.parents[1]/"cavity_triplet"/"online"/"dock.py"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay(work, engine, shape, pro, dest):
    # Pin both the controller and its children, including simultaneous imports.
    core = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {core})
    configs = json.loads((work/"inputs.json").read_text())["cells"]
    dest.mkdir(exist_ok=True,parents=True)
    rows = []
    for j, src in enumerate(sorted((work/"banks").glob("*.source.npz"))):
        stem = src.name.removesuffix(".source.npz")
        name, seed = stem.rsplit("_",1)
        cell = configs[name]
        arms = ["pro","shape","true","perm"]
        for arm in arms[j%4:]+arms[:j%4]:
            tag = f"{stem}_{arm}"
            out = dest/(tag+".pdbqt")
            runner, exe, mode = ((SHAPE_ADAPTER, shape, "matched") if arm=="shape"
                                 else (HERE/"live.py",engine,arm))
            if arm == "pro":
                exe = pro
            prefix = [pro] if arm=="pro" else [sys.executable, runner, "--engine", exe, "--arm", mode, "--"]
            cmd = [*prefix, "--ligand", work/cell["ligand"], "--ligand_bonds", work/cell["bonds"],
                   "--receptor", work/cell["receptor"], "--cpu", 1, "--seed", seed,
                   "--max_evals", 25000, "--exhaustiveness", 2, "--num_modes", 9,
                   "--energy_range", 100, "--pose_selection", "late", "--search_counters", "--out", out]
            for axis, center in zip("xyz",cell["center"],strict=True):
                cmd.extend(["--center_"+axis, center, "--size_"+axis, cell["size"]])
            log = out.with_suffix(".log")
            with log.open("w") as f:
                before = resource.getrusage(resource.RUSAGE_CHILDREN)
                t0 = time.perf_counter()
                subprocess.run(list(map(str,cmd)), env=os.environ | {"OPENBLAS_NUM_THREADS":"1","OMP_NUM_THREADS":"1"},
                               stdout=f, stderr=subprocess.STDOUT, check=True)
                wall_s = time.perf_counter()-t0
                after = resource.getrusage(resource.RUSAGE_CHILDREN)
            raw = work/"raw"
            if (raw/(tag+".pdbqt")).read_bytes() != out.read_bytes():
                raise ValueError(f"{tag}: published pool differs")
            def counters(path):
                return sorted(line for line in path.read_text().splitlines()
                              if line.startswith(("Search trajectory:", "Search counters:", "Hunt chain:")))
            if counters(raw/(tag+".log")) != counters(log):
                raise ValueError(f"{tag}: local solver counts differ")
            if arm in ("true","perm"):
                def winner(path):
                    return next(l for l in path.read_text().splitlines() if l.startswith("field_rank "))
                assert winner(raw/(tag+".inject"))==winner(log)
            rows.append(dict(case=stem, arm=arm, wall_s=wall_s,
                             cpu_s=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                             pool_sha256=digest(out), engine_sha256=digest(exe), pinned_cpu=core))
    with (dest/"cost.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,rows[0]); w.writeheader(); w.writerows(rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("pose_work","rank_engine","shape_engine","pro_engine","dest"):
        p.add_argument(name,type=lambda arg:Path(arg).resolve())
    a=p.parse_args()
    replay(a.pose_work,a.rank_engine,a.shape_engine,a.pro_engine,a.dest)
