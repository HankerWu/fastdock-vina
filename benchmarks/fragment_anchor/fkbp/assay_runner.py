"""Prepared-input single-pass AD4 comparisons. Only check_panel.py uses Ki and crystal data."""
import argparse
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import sys
import tarfile
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SERIES = ROOT / "benchmarks/cubic_grid"
RECEPTORS = ("1FKG", "1FKH")
ARMS = ("official", "base", "whole", "fragment", "sham", "high", "official_high")
LIVE = ("whole", "fragment", "sham")
SEED = 1047293
OFFICIAL_SHA = "f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def prep(work):
    sys.path.insert(0, str(ROOT / "benchmarks/receptor_linkage/aurora/internal_loops/rigid_fkbp"))
    from run_panel import prepare
    hashes = prepare(work)
    archive = SERIES / "human_panel/prepared.tar.xz"
    hashes[str(archive.relative_to(ROOT))] = sha(archive)
    with tarfile.open(archive) as t:
        for rec in RECEPTORS:
            t.extract(rec + "/receptor.pdbqt", work, filter="data")
    from rdkit import Chem
    for ligand in sorted((work / "affinity").glob("*.pdbqt")):
        lines = ligand.read_text().splitlines()
        mol = Chem.MolFromSmiles(next(s.split()[2] for s in lines if s.startswith("REMARK SMILES ") and not s.startswith("REMARK SMILES IDX")))
        mapped = {}
        for line in lines:
            if line.startswith("REMARK SMILES IDX "):
                vals = list(map(int, line.split()[3:]))
                mapped.update(zip(vals[::2], vals[1::2], strict=True))
        assert set(mapped) == set(range(1, 1+mol.GetNumAtoms()))
        ligand.with_suffix(".bonds").write_text("".join(
            f"{mapped[b.GetBeginAtomIdx()+1]} {mapped[b.GetEndAtomIdx()+1]}\n" for b in mol.GetBonds()))
    paths = sorted([*work.glob("affinity/*"), *work.glob("*/receptor.pdbqt"), *work.glob("maps/*/*")])
    (work / "inputs.json").write_text(json.dumps(dict(archives=hashes,files={str(p.relative_to(work)):sha(p) for p in paths}), indent=2)+"\n")


def live_one(work, rec, n, arm, exe, stem, conf):
    import numpy as np
    from scipy.spatial.transform import Rotation
    sys.path[:0] = [str(ROOT / "benchmarks/cavity_field_rank/engine"),str(HERE.parent / "engine"),
                    str(ROOT / "benchmarks/cavity_field_rank/fkbp"),str(ROOT/"tools")]
    from anchors import place
    from live import open_bank
    from frozen_assay import atom_types, field_score, old
    from pocket_triangle import Pocket
    pipe = stem.with_suffix(".fifo")
    os.mkfifo(pipe, mode=0o600)
    qfile, xfile = stem.with_suffix(".q"), stem.with_suffix(".source")
    event = stem.with_suffix(".event")
    event.write_text("1 1 0 0 0 1 0 0 0\n")
    env = dict(os.environ, FD_SOURCE_TRAJECTORY="1", FD_SOURCE_STEP="1",
               FD_SOURCE_Q=str(qfile),FD_SOURCE_XYZ=str(xfile), FD_RIGID_INTERVENTION=str(event),
               FD_RIGID_LOG=str(stem.with_suffix(".stages")), FD_RIGID_COORDS=str(stem.with_suffix(".xyz")),
               FD_RIGID_HISTORY=str(stem.with_suffix(".history")),FD_FRAGMENT_TRANSFORM=str(pipe))
    cmd = [str(exe), "--config", str(conf.relative_to(work)), "--search_counters"]
    try:
        with subprocess.Popen(cmd, cwd=work, env=env) as child:
            try:
                with open_bank(pipe, child) as stream:
                    header, line = xfile.read_text().splitlines()
                    src_serial = np.array(header[1:].split(),int)
                    x0 = np.array(line.split()[3:],float).reshape(-1,3)
                    order = np.argsort(src_serial)
                    serial, x = src_serial[order], x0[order]
                    q = np.loadtxt(qfile)[2:]
                    ligand = work / "affinity" / f"{n}.pdbqt"
                    with (work/"maps"/rec/"receptor.C.map").open() as f:
                        header = [next(f).split() for _ in range(6)]
                    center, size = np.array(header[5][1:],float), np.array(header[4][1:],int)*float(header[3][1])
                    pocket = Pocket.from_pdbqt(work/rec/"receptor.pdbqt",center,size)
                    bank = place(pocket,x,serial,ligand,ligand.with_suffix(".bonds"),SEED,arm)
                    rot = bank["matched_rotation"]
                    poses = (x-x.mean(0)) @ rot.transpose(0,2,1) + bank["center"][:,None]
                    real, _ = atom_types(ligand, serial)
                    fields = {t:old.read_map(work/"maps"/rec/("receptor."+t+".map")) for t in sorted(set(real))}
                    scores, outside = field_score(poses,real,fields)
                    assert len(scores) == 1024
                    i = int(np.argmin(scores))
                    r = Rotation.from_matrix(rot[i])
                    orig = r.apply(q[:3]-x.mean(0)) + bank["center"][i]
                    quat = (r*Rotation.from_quat(q[[4,5,6,3]])).as_quat()[[3,0,1,2]]
                    transform = " ".join(format(v,".17g") for v in np.r_[orig,quat])+"\n"
                    stream.write(transform)
                    stream.flush()
                # Compact audit data; transforms can also be reconstructed from ordinary inputs and .source.
                (stem.with_suffix(".transform")).write_text(transform)
                np.savez(stem.with_suffix(".npz"),site_triangles=bank["site_triangles"],anchor=bank["ligand_triangle"],anchor_serials=bank["anchor_serials"],
                         component_serials=bank["component_serials"],fallback=bank["anchor_fallback"],
                         true_score=scores,center=bank["center"],rotation=rot,shape=bank["matched_score"])
                with stem.with_suffix(".rank").open("w") as f:
                    json.dump(dict(winner=i,queries=int(poses.shape[0]*poses.shape[1]),outside=int(outside.sum()),
                              pocket_grid_points=pocket.n_grid,pocket_burial_points=pocket.n_burial,
                              shape_queries=int(poses.shape[0]*poses.shape[1]), ordered_triangles=len(pocket.triangles)),f)
                status = child.wait()
                if status: raise subprocess.CalledProcessError(status,cmd)
            except BaseException:
                if child.poll() is None: child.terminate()
                raise
    finally:
        pipe.unlink(missing_ok=True)


def time_one(cmd, work, path):
    env = {k:v for k,v in os.environ.items() if not k.startswith("FD_")}
    def cpu():
        u=resource.getrusage(resource.RUSAGE_CHILDREN)
        return u.ru_utime+u.ru_stime
    start, before = time.perf_counter(), cpu()
    with path.open("w") as f:
        subprocess.run(cmd, cwd=work, env=env, stdout=f,stderr=subprocess.STDOUT,check=True)
    return time.perf_counter()-start,cpu()-before


def run(work, official, engine, counter, only):
    assert sha(official) == OFFICIAL_SHA
    for j,(rec,n) in enumerate(itertools.product(RECEPTORS,range(4,12))):
        for arm in ARMS[j%len(ARMS):]+ARMS[:j%len(ARMS)]:
            if arm not in only: continue
            stem = work/"pools"/arm/f"{rec}_{n}"
            stem.parent.mkdir(parents=True,exist_ok=True)
            conf, out = stem.with_suffix(".conf"),stem.with_suffix(".pdbqt")
            config = dict(ligand=f"affinity/{n}.pdbqt",maps=f"maps/{rec}/receptor",scoring="ad4",num_modes=9,
                          energy_range=100,cpu=1,exhaustiveness=3,max_evals=48000 if "high" in arm else 24000,
                          seed=SEED,out=str(out.relative_to(work)))
            conf.write_text("".join(f"{k} = {v}\n" for k,v in config.items()))
            is_off = arm.startswith("official")
            program = official if is_off else engine
            cmd = [str(program),"--config",str(conf.relative_to(work))]
            if arm in LIVE:
                cmd = [sys.executable,str(Path(__file__).resolve()),"live",str(work),rec,str(n),arm,str(engine)]
            elif not is_off:
                cmd.append("--search_counters")
            log = stem.with_suffix(".log")
            wall,cpu = time_one(cmd,work,log)
            replay_wall = replay_cpu = 0.0
            if is_off:
                log = stem.with_suffix(".replay.log")
                alt = stem.with_suffix(".replay")
                replay_wall,replay_cpu = time_one([str(counter),"--config",str(conf.relative_to(work)),"--out",
                                                 str(alt.relative_to(work))],work,log)
                assert alt.read_bytes() == out.read_bytes()
                alt.unlink()
            output = log.read_text()
            records = re.findall(r"Search trajectory:[^\n]+",output)
            assert len(records)==3
            def count(key):
                return sum(int(re.search(r"\b"+key+r"=(\d+)",r)[1]) for r in records)
            if is_off:
                channels = sum(int(v) for v in re.findall(r"GRID_CHANNEL derivative=\d+ n=(\d+)",output))
            else:
                channels = 3*sum(int(v) for v in re.findall(r"FIELD_COUNT combined=1 derivative=\d+ positions=(\d+)",output))
            assert channels>0
            ranks = json.loads(stem.with_suffix(".rank").read_text()) if arm in LIVE else {}
            row = dict(arm=arm,receiver=rec,compound=n,seed=SEED,modes=out.read_text().count("ENDMDL"),
                sha256=sha(out),executable_sha256=sha(program),counter_sha256=sha(counter) if is_off else "",
                hunt_evals=count("hunt_evals"),second_evals=count("second_evals"),steps=count("hunt_calls"),
                coarse_channels=ranks.get("queries",0),downstream_channels=channels,wall_s=wall,cpu_s=cpu,
                official_counter_replay_wall_s=replay_wall,official_counter_replay_cpu_s=replay_cpu)
            stem.with_suffix(".timing").write_text(json.dumps(row)+"\n")
            print(arm,rec,n,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="action",required=True)
    sub.add_parser("prepare").add_argument("work",type=Path)
    runp=sub.add_parser("run")
    for arg in ("work","official","engine","counter"):
        runp.add_argument(arg,type=Path)
    runp.add_argument("--arms",nargs="+",choices=ARMS,default=ARMS)
    lp=sub.add_parser("live")
    for arg,t in (("work",Path),("rec",str),("compound",int),("arm",str),("exe",Path)): lp.add_argument(arg,type=t)
    a=p.parse_args()
    if a.action=="prepare": prep(a.work.resolve())
    elif a.action=="run": run(a.work.resolve(),a.official.resolve(),a.engine.resolve(),a.counter.resolve(),a.arms)
    else:
        stem=a.work/"pools"/a.arm/f"{a.rec}_{a.compound}"
        live_one(a.work,a.rec,a.compound,a.arm,a.exe,stem,stem.with_suffix(".conf"))
