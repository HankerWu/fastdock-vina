"""Actual single-event rigid forks and unmodified search controls, one CPU per process."""
import argparse
import csv
import hashlib
import json
import os
import re
import resource
import subprocess
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "receptor_linkage/aurora/internal_loops"
SEEDS = (91813, 91823)
COUNT_KEYS = ("hunt_calls", "hunt_evals", "second_calls", "second_evals", "explicit_calls",
              "explicit_evals", "late_prescore")


def digest(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


def table(filename, rows):
    with filename.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def systems(p):
    # Reuse the published, source-graph-checked inputs; no query reference or
    # evaluation coordinates are loaded anywhere in the blind docking driver.
    cells = json.loads((OLD / "rigid_coupled/results/inputs.json").read_text())["inputs"]
    for name, c in cells.items():
        if name == "82I":
            paths = [p/"brd/6YQW/ligand.pdbqt", OLD/"rigid_coupled/results/82I.bonds",
                     p/"brd/6Y7L/receptor.pdbqt"]
        else:
            edges = "at9283_chemical.bonds" if name == "at9283" else "danusertib_native.bonds"
            paths = [p/f"aurora/prepared/{name}.pdbqt", p/"aurora/raw"/edges,
                     p/"aurora/template1OL7/trans_model/released/receptor.pdbqt"]
        for key, path in zip(("ligand", "bonds", "receptor"), paths, strict=True):
            if digest(path) != c[key]:
                raise ValueError(f"{name} {key} identity")
            c[key] = path
    return cells


def invoke(p, costs, stem, binary, c, seed, cap=25000, official=False, env=None):
    pool = p/"raw"/(stem+".pdbqt")
    log = pool.with_suffix(".log")
    args = [binary, "--ligand", c["ligand"], "--receptor", c["receptor"]]
    for axis, center in zip("xyz", c["center"], strict=True):
        args.extend([f"--center_{axis}", center, f"--size_{axis}", c["size"]])
    args.extend(["--seed", seed, "--exhaustiveness", 2, "--cpu", 1, "--max_evals", cap,
                 "--num_modes", 9, "--energy_range", 100, "--out", pool])
    if not official:
        args.extend(["--ligand_bonds", c["bonds"], "--pose_selection", "late", "--search_counters"])
    start = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    clean = {k: v for k, v in os.environ.items() if not k.startswith("FD_")}
    with log.open("w") as f:
        subprocess.run(list(map(str,args)), env=clean | {k: str(v) for k,v in (env or {}).items()},
                       stdout=f, stderr=subprocess.STDOUT, check=True)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    costs.append(dict(stem=stem, wall_s=time.perf_counter()-start,
                      cpu_s=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                      pool_sha256=digest(pool), **counts(log)))
    return pool


def counts(log):
    s = log.read_text()
    return {k: sum(map(int,re.findall(r"\b"+k+r"=(\d+)",s))) for k in COUNT_KEYS}


def source_record(p, name, c, seed, tr, step, xyz, q, serials):
    source = f"{name}_{seed}_t{tr}_s{step}"
    file = p/"sources"/(source+".npz")
    np.savez_compressed(file, xyz=xyz, q=q, serials=serials)
    return dict(source=source,drug=name,seed=seed,trajectory=tr,step=step,
                receptor=c["receptor"],ligand=c["ligand"],bonds=c["bonds"],
                **dict(zip(("center_x","center_y","center_z"),c["center"],strict=True)),
                size=c["size"],npz=file)


def sources(p, cells, binary):
    data=np.load(OLD/"proposal_crossover/results/complete_on_complete.npz")
    rows=[]
    for i,(tr,step) in enumerate(zip(data["trajectory"],data["step"],strict=True)):
        if step == 1 or step % 100 == 0:
            rows.append(source_record(p,"at9283",cells["at9283"],61231,tr,step,
                                      data["start_xyz"][i],data["generalized"][i],data["serials"]))
    costs=[]
    for name, c in cells.items():
        for seed in SEEDS:
            stem=f"{name}_{seed}_source"
            xyz_path=p/"raw"/(stem+".xyz")
            q_path=xyz_path.with_suffix(".q")
            invoke(p,costs,stem,binary,c,seed,env=dict(FD_SOURCE_STEP=1,FD_SOURCE_TRAJECTORY=1,
                   FD_SOURCE_XYZ=xyz_path,FD_SOURCE_Q=q_path))
            header, line=xyz_path.read_text().splitlines()
            serials=np.array(header.split()[1:],int)
            xyz=np.array(line.split()[3:],float).reshape(-1,3)
            q=np.loadtxt(q_path)
            assert tuple(q[:2]) == (1,1) and line.split()[:3] == ["source","1","1"]
            rows.append(source_record(p,name,c,seed,1,1,xyz,q[2:],serials))
    assert len(rows)==24
    table(p/"sources.csv",rows)
    table(p/"source_cost.csv",costs)


def controls(p,cells,a):
    expected=json.loads((OLD/"rigid_coupled/results/inputs.json").read_text())["binary_sha256"]
    costs=[]
    checks=[]
    if digest(a.official) != expected["official"]:
        raise ValueError("Expected the official release executable")
    # Compiler paths can alter the counter binary. Its modes must match the
    # pinned release byte-for-byte on every invocation; report its actual hash.
    for name,c in cells.items():
        for seed in SEEDS:
            for arm,cap in (("pro",25000),("pro_double",50000),("official",25000),("official_double",50000)):
                stem=f"{name}_{seed}_{arm}"
                official=arm.startswith("official")
                pool=invoke(p,costs,stem,a.official if official else a.pro,c,seed,cap,official)
                if arm=="pro":
                    assert pool.read_bytes()==(p/"raw"/f"{name}_{seed}_source.pdbqt").read_bytes()
                if official:
                    check=invoke(p,checks,stem+"_counter",a.counter,c,seed,cap,True)
                    assert check.read_bytes()==pool.read_bytes() and checks[-1]["hunt_evals"]>0
                    costs[-1].update({k:checks[-1][k] for k in COUNT_KEYS})
    table(p/"control_cost.csv",costs)
    table(p/"counter_cost.csv",checks)
    inputs={name:{k:(digest(c[k]) if k in ("ligand","receptor","bonds") else c[k])
                  for k in ("ligand","receptor","bonds","center","size")} for name,c in cells.items()}
    (p/"input_hashes.json").write_text(json.dumps(dict(systems=inputs,binaries={
        key:digest(getattr(a,key)) for key in ("official","counter","pro","observer")}),indent=2)+"\n")


def branches(p,cells,binary, timed=False, posthoc=False):
    sources={r["source"]:r for r in csv.DictReader((p/"sources.csv").open())}
    costs=[]
    selections=list(csv.DictReader((p/("posthoc_selections.csv" if posthoc else "selections.csv")).open()))
    assert len(selections)==(4 if posthoc else 2*len(sources))
    for row in selections:
        source=row["source"]
        r=sources[source]
        stem=source+"_"+row["arm"]
        s=p/"raw"/(stem+".input")
        stage=s.with_suffix(".stages")
        log=s.with_suffix(".inject")
        history=s.with_suffix(".history")
        q=[float(row[k]) for k in ("qx","qy","qz","qw","qi","qj","qk")]
        s.write_text(" ".join([r["trajectory"],r["step"],*[f"{v:.17g}" for v in q]])+"\n")
        env = dict(FD_RIGID_INTERVENTION=s, FD_RIGID_LOG=log)
        run_stem=stem
        if timed:
            run_stem += "_timed"
            log=s.with_suffix(".timed_inject")
            env["FD_RIGID_LOG"]=log
        else:
            env.update(FD_RIGID_COORDS=stage, FD_RIGID_HISTORY=history)
        pool=invoke(p,costs,run_stem,binary,cells[r["drug"]],int(r["seed"]),env=env)
        if timed:
            assert pool.read_bytes()==(p/"raw"/(stem+".pdbqt")).read_bytes()
        original=np.array(log.read_text().splitlines()[0].split()[1:],float)
        np.testing.assert_array_equal(original,np.load(r["npz"])["q"])
    table(p/("posthoc_cost.csv" if posthoc else ("prediction_cost.csv" if timed else "branch_cost.csv")),costs)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=["sources","controls","branches","timing","posthoc"])
    parser.add_argument("work",type=lambda v:Path(v).resolve())
    for k in ("observer","official","counter","pro"):
        parser.add_argument("--"+k,type=lambda v:Path(v).resolve())
    a=parser.parse_args()
    p=a.work
    (p/"raw").mkdir(exist_ok=True)
    (p/"sources").mkdir(exist_ok=True)
    c=systems(p)
    if a.action=="sources": sources(p,c,a.observer)
    if a.action=="controls": controls(p,c,a)
    if a.action in ("branches","timing","posthoc"):
        branches(p,c,a.observer,a.action=="timing",a.action=="posthoc")
