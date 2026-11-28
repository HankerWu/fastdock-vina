"""Cross the two source-pinned hunt solvers on frozen historical proposals and verify native endpoints."""
import argparse
import csv
import hashlib
import json
import lzma
import os
from pathlib import Path
import sys
import tarfile
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/"stage_local_structure"))
from external import structures, run_panel
sys.path.insert(0, str(HERE.parent/"bfgs_source_replay"))
from analyze import native
from archive import restore

ARMS = ("stable", "candidate")
TRACE = HERE.parent/"stage_local_hunt/coverage"
TRIAL = ("index iteration attempt alpha slope fd4 fd5 initial final predicted actual "
         "translation rotation torsion rms crossed outside0 outside1 accepted").split()


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def records(path, data=None):
    if data is not None:
        path.write_bytes(lzma.compress(data, preset=3))
    return [line.split("\t") for line in lzma.decompress(path.read_bytes()).decode().splitlines()]


def csvfile(path, rows, live):
    if live:
        run_panel.csvfile(path, rows)
    else:
        with path.open() as stream:
            assert list(csv.DictReader(stream)) == [{k: str(v) for k, v in r.items()} for r in rows], path


def fold(trajectories):
    return [step for traj in trajectories for step in traj]


def forks(runs, case, table):
    probes = {arm: [] for arm in ARMS}
    offsets = [0, 0]
    for trajectory, pair in enumerate(zip(*(runs[arm] for arm in ARMS)), 1):
        def first_difference(key):
            for i, (a, b) in enumerate(zip(*pair)):
                if key == "A":
                    left, right = bool(a[key]), bool(b[key])
                else:
                    fields = slice(3,4) if key == "P" else slice(2,4)
                    left, right = a[key][0][fields], b[key][0][fields]
                if left != right:
                    return i
            return None
        split = {k: first_difference(k) for k in ("S", "P", "H", "A")}
        i = split["H"]
        assert i is not None
        row = dict(case=case, trajectory=trajectory,
                   **{"first_"+k: "" if v is None else v for k,v in split.items()})
        for j, arm in enumerate(ARMS):
            t = pair[j]
            idx = offsets[j]+i
            probes[arm].append(idx)
            offsets[j] += len(t)
            row.update({arm+"_index": idx, arm+"_source":t[i]["S"][0][3],
                        arm+"_proposal": t[i]["P"][0][3], arm+"_hunt": t[i]["H"][0][3],
                        arm+"_energy":t[i]["H"][0][2], arm+"_accepted": bool(t[i]["A"])})
        table.append(row)
    return probes


def argv(binary, inputs, q, r, seed, target):
    cmd = [str(binary), "--ligand", str(inputs/q["pdb"]/"ligand.pdbqt"), "--receptor",
           str(inputs/r["pdb"]/"receptor.pdbqt"), "--seed", str(seed), "--cpu", "1",
           "--exhaustiveness", "3", "--max_evals", "12000", "--num_modes", "9", "--energy_range", "100",
           "--spacing", "0.375", "--out", str(target)]
    for axis, center in zip("xyz", r["center"]):
        cmd += ["--size_"+axis, "24", "--center_"+axis, str(center)]
    return cmd


def analyze_trials(key, data, selected, result):
    from itertools import groupby
    values = [dict(zip(TRIAL, map(float, row.split(",")))) for row in data.decode().splitlines()]
    assert sorted({int(r["index"]) for r in values}) == sorted(selected)
    groups = []
    for ((idx, direction), group) in groupby(values, key=lambda r: (int(r["index"]), int(r["iteration"]))):
        trials = list(group)
        assert all(r["attempt"] == j and r["alpha"] == 0.5**j for j,r in enumerate(trials))
        assert all(not r["accepted"] for r in trials[:-1])
        last = trials[-1]
        assert last["accepted"] or len(trials) == 10
        assert last["slope"] < 0 and all(np.isfinite(last[k]) for k in TRIAL)
        if last["accepted"]:
            assert last["actual"] < 0
        groups.append((idx, trials[0], last))
    for idx, g in groupby(groups, key=lambda item: item[0]):
        path = list(g)
        calls = 1+sum(int(last["attempt"])+1 for _, _, last in path)
        assert calls == int(result[idx][4]), (key, idx, calls, result[idx][4])
        assert sum(not last["accepted"] for _,_,last in path) == int(result[idx][5])
    rows = []
    for passed in (True,False):
        g = [v for v in groups if bool(v[2]["accepted"])==passed]
        if not g:
            continue
        terminal = [last for _,_,last in g]
        b = sum(r["fd4"]<0 and r["fd5"]<0 and (r["crossed"]>0 or r["rms"]>=.375) for r in terminal)
        row = dict(key=key, end="pass" if passed else "exhausted", directions=len(g),
                   trials=sum(int(r["attempt"]+1) for r in terminal), diagnostic_Egrad=4*len(g),
                   negative_both=sum(r["fd4"]<0 and r["fd5"]<0 for r in terminal),
                   changed_cell=sum(r["crossed"]>0 for r in terminal),
                   outside=sum(r["outside1"]>0 for r in terminal),
                   registered_B=b)
        for name in ("rms","translation","rotation","torsion","actual","predicted"):
            row["first_"+name+"_median"] = float(np.median([first[name] for _,first,_ in g]))
            row["last_"+name+"_median"] = float(np.median([r[name] for r in terminal]))
        err = [abs(r["fd5"]-r["slope"])/max(1,abs(r["slope"]),abs(r["fd5"])) for r in terminal]
        row["fd5_scaled_error_p50"],row["fd5_scaled_error_p95"] = map(float,np.quantile(err,[.5,.95]))
        rows.append(row)
    return rows


def main(work,out,builds,free_builds,historical):
    live=builds is not None
    if live:
        out.mkdir(parents=True,exist_ok=False)
    if not historical:
        historical=work/"history"
        restore(HERE.parent/"bfgs_source_replay/raw",historical)
    ext=work/"external"
    with tarfile.open(HERE.parent/"stage_local_structure/external/inputs.tar.xz") as archive:
        archive.extractall(ext,filter="data")
    e, omitted=structures(ext)
    assert not omitted
    histfix=historical/"fixtures/reciprocal"
    f=[json.loads((histfix/c/"preparation.json").read_text()) for c in ("2NNQ","3FR2")]
    costs=[]
    if not live:
        old=json.loads((out/"costs.json").read_text())
    fork_rows, bank_rows, trial_rows=[],[],[]
    for inputs,q,r,seed,drivers in [(histfix,*f,85237,("e7031e9","2326430")),
                                     (ext,*e,85231,ARMS),(ext,*e[::-1],85231,ARMS)]:
        case=f"{q['pdb']}-{r['pdb']}-{seed}"
        hist=historical/case
        trajectories={}
        if inputs==ext:
            for arm in ARMS:
                raw=lzma.decompress((TRACE/f"reciprocal-{case}-{arm}.xz").read_bytes())
                nativefile=work/f"{case}-{arm}.native"
                if live and free_builds:
                    binary=free_builds/arm/"build/vina"
                    pose=work/f"{case}-{arm}.pdbqt"
                    cmd=argv(binary,inputs,q,r,seed,pose)
                    os.environ.pop("FD_BFGS_SHADOW",None)
                    os.environ["FD_BFGS_TRACE"]=str(nativefile)
                    stats=run_panel.execute(cmd,out/f"{case}-{arm}.free.log",180)
                    assert stats["status"]=="ok" and nativefile.read_bytes()==raw
                    assert pose.read_bytes()==(TRACE/f"reciprocal-{case}-{arm}.pdbqt").read_bytes()
                    costs.append(dict(key=f"{case}-{arm}-free",binary=run_panel.digest(binary),
                                      trace=digest(raw),argv=cmd,**stats))
                else:
                    nativefile.write_bytes(raw)
                trajectories[arm]=native(nativefile)[0]
            extra=forks(trajectories,case,fork_rows)
        else:
            trajectories={arm:native(hist/f"{arm}.native")[0] for arm in drivers}
            extra={arm:[] for arm in drivers}
        ref=run_panel.replay.reciprocal.reference(q,r)
        serials=[int(l[6:11]) for l in (inputs/q["pdb"]/"ligand.pdbqt").read_text().splitlines()
                 if l.startswith(("ATOM  ","HETATM")) and l[77:].strip() not in {"H","HD","HS"}]
        order=[serials.index(i) for i in q["serial"]]
        def rms(row):
            xyz=np.fromstring(row[8],dtype=int,sep=" ").reshape(-1,3)[order]/1000
            return float(np.sqrt(np.square(xyz[None]-ref).sum(-1).mean(-1).min()))
        for driver in drivers:
            free=fold(trajectories[driver])
            source=work/f"{case}-{driver}.inputs"
            source.write_text("".join(f"{s['P'][0][3]}\n" for s in free))
            selected=set(range(0,len(free),97))|set(extra[driver])
            index=source.with_suffix(".indices")
            index.write_text("".join(f"{i}\n" for i in sorted(selected)))
            if hist.exists():
                assert source.read_bytes()==(hist/f"{driver}.inputs").read_bytes()
            for arm in ARMS:
                key=f"{case}-{driver}-{arm}"
                trace=work/f"{key}.shadow"
                trials=work/f"{key}.trials"
                archive=out/f"{key}.shadow.xz"
                packed=out/f"{key}.trials.xz"
                if live:
                    binary=builds/arm/"build/vina"
                    cmd=argv(binary,inputs,q,r,seed,trace.with_suffix(".unused"))
                    os.environ.update(FD_BFGS_SHADOW=str(source),FD_BFGS_TRACE=str(trace),
                                      FD_FORK_INDICES=str(index),FD_FORK_TRIALS=str(trials))
                    stats=run_panel.execute(cmd,out/f"{key}.log",180)
                    if stats["status"]!="ok":
                        raise ValueError((key,stats))
                    result=records(archive,trace.read_bytes())
                    packed.write_bytes(lzma.compress(trials.read_bytes(),preset=3))
                    costs.append(dict(key=key,argv=cmd,binary=run_panel.digest(binary),
                            inputs={n:run_panel.digest(path) for n,path in {
                                "states":source,"indices":index,"receptor":inputs/r["pdb"]/"receptor.pdbqt",
                                "ligand":inputs/q["pdb"]/"ligand.pdbqt"}.items()},
                            shadow_sha256=run_panel.digest(archive),trials_sha256=run_panel.digest(packed),**stats))
                else:
                    receipt=next(v for v in old if v["key"]==key)
                    assert run_panel.digest(source)==receipt["inputs"]["states"] and run_panel.digest(index)==receipt["inputs"]["indices"]
                    assert run_panel.digest(archive)==receipt["shadow_sha256"] and run_panel.digest(packed)==receipt["trials_sha256"]
                    result=records(archive)
                assert len(result)==len(free) and all(row[1]==str(i) for i,row in enumerate(result))
                if hist.exists():
                    expect=([line.split("\t") for line in (hist/f"{driver}-ed0d42d.shadow").read_text().splitlines()]
                        if arm=="stable" else records(HERE.parent/"stage_local_hunt/locked"/f"{case}-{driver}.shadow.xz"))
                    assert all(a==b for a,b in zip(result,expect))
                elif driver==arm:
                    assert all(a[2:7]==s["H"][0][2:7] and a[8]==s["H"][0][8] for a,s in zip(result,free))
                trial_rows.extend(analyze_trials(key,lzma.decompress(packed.read_bytes()),selected,result))
                errors=[rms(row) for row in result]
                bank_rows.append(dict(case=case,driver=driver,solver=arm,proposals=len(result),
                     calls=sum(int(row[4]) for row in result), exhausted=sum(int(row[5]) for row in result),
                     selected=len(selected),near=sum(x<=2 for x in errors),best_A=min(errors),
                     near_on_driver_accepted=sum(x<=2 and bool(s["A"]) for x,s in zip(errors,free))))
                print(bank_rows[-1],flush=True)
    for name,rows in (("forks",fork_rows),("banks",bank_rows),("directions",trial_rows)):
        csvfile(out/f"{name}.csv",rows,live)
    if live:
        (out/"costs.json").write_text(json.dumps(costs,indent=2)+"\n")
    print("Verified archived same-solver outcomes, paired first forks and all selected directions.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    p.add_argument("--builds",type=Path)
    p.add_argument("--free-builds",type=Path)
    p.add_argument("--history",type=Path)
    a=p.parse_args()
    with tempfile.TemporaryDirectory() as work:
        main(Path(work),a.out.resolve(),a.builds,a.free_builds,a.history)
