"""Independently verify candidate tangents and the archived outputs for the six locked libraries."""
import argparse
import csv
import hashlib
import json
import lzma
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import numpy as np
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/"bfgs_source_replay"))
from instrument_candidate import edit,instrument
from verify import reference
from run_replay import run_panel
CANDIDATE="0bcae6ef09857f97656291bb4821d20b3b619bf5"
def sha(data):
    return hashlib.sha256(data).hexdigest()
def prepare(path):
    path.mkdir()
    stage=subprocess.run(["git","-C",REPO,"archive",CANDIDATE,"src","tools/build.py"],stdout=subprocess.PIPE,check=True)
    import io
    with tarfile.open(fileobj=io.BytesIO(stage.stdout)) as tar:
        tar.extractall(path,filter="data")
    instrument(path)
    lib=path/"src/lib"
    edit(lib/"bfgs_trace.h","local(m, p, ig, pose, cap, cost);","local(m, p, ig, pose, cap, cost, nullptr, 0, true);")
    shutil.copyfile(HERE/"audit_direction.h",lib/"audit_direction.h")
    edit(lib/"bfgs.h","struct bfgs_result {",'#include "audit_direction.h"\nstruct bfgs_result {')
    edit(lib/"bfgs.h","if (hunt_direction) limit_hunt_direction(x, gradient, direction);", """if (hunt_direction) {
                const flv raw = direction;
                limit_hunt_direction(x, gradient, direction);
                static direction_audit diagnostic;
                diagnostic.check(x, gradient, raw, direction);
            }""")
def evaluate(history,work,binary):
    output=HERE.parent/"stage_local_hunt/locked"
    original=list(csv.DictReader((output/"locked.csv").open()))
    receipts=json.loads((output/"costs.json").read_text())
    result=[]
    for item in receipts:
        case,driver=item["case"],item["driver"]
        q,r,seed=case.split("-")
        stem=f"{case}-{driver}"
        old=history/case
        fixture=history/"fixtures/reciprocal"
        query,receiver=(json.loads((fixture/key/"preparation.json").read_text()) for key in (q,r))
        ref=reference(query,receiver)
        serials=[int(l[6:11]) for l in (fixture/q/"ligand.pdbqt").read_text().splitlines()
                 if l.startswith(("ATOM","HETATM")) and int(l[6:11]) in query["serial"]]
        order=[serials.index(i) for i in query["serial"]]
        def rmsd(x):
            xyz=np.fromstring(x[8],sep=" ").reshape(-1,3)[order]*0.001
            delta=xyz-ref
            return np.sqrt(np.einsum("sij,sij->s",delta,delta).min()/len(order))
        states=old/(driver+".inputs")
        assert sha(states.read_bytes())==item["states_sha256"]
        cmd=[str(binary),"--ligand",str(fixture/q/"ligand.pdbqt"),"--receptor",str(fixture/r/"receptor.pdbqt"),
             "--cpu","1","--seed",seed,"--max_evals","12000","--exhaustiveness","3","--num_modes","9",
             "--energy_range","100","--spacing","0.375","--out",str(work/(stem+".unused"))]
        for axis,center in zip("xyz",receiver["center"]):
            cmd+=["--size_"+axis,"24","--center_"+axis,str(center)]
        audit=work/(stem+".audit")
        trace=work/(stem+".shadow")
        os.environ.update(FD_AUDIT_DIRECTION=str(audit),FD_BFGS_SHADOW=str(states),FD_BFGS_TRACE=str(trace))
        stats=run_panel.execute(cmd,work/(stem+".log"),180)
        assert stats["status"]=="ok"
        data=trace.read_bytes()
        assert sha(data)==item["raw_sha256"]
        assert data==lzma.decompress((output/(stem+".shadow.xz")).read_bytes())
        events=[line.split("\t") for line in data.decode().splitlines()]
        assert all(np.isfinite(float(e[10])) and float(e[10])<0 and int(e[6])==0 for e in events)
        for i,e in enumerate(events):
            assert e[0]=="H" and int(e[1])==i
            if int(e[5]):
                assert np.isfinite(float(e[11])) and float(e[11])<0 and 0<=float(e[12])<=np.pi/512 +1e-12
        for label,rev in (("stable","ed0d42d"),("common_scalar","0bb827e")):
            coords=[row.split("\t") for row in (old/(driver+"-"+rev+".shadow")).read_text().splitlines()]
            assert len(coords)==len(events)==len(states.read_text().splitlines())
            assert all(a[9]==b[9] for a,b in zip(coords,events))
            expect=next(x for x in original if x["case"]==case and x["driver"]==driver and x["control"]==label)
            values=np.array([rmsd(f) for f in events]); historical=np.array([rmsd(f) for f in coords])
            assert abs(values.min()-float(expect["candidate_best_A"]))<1e-9
            assert abs(historical.min()-float(expect["control_best_A"]))<1e-9
            assert sum(int(x[4]) for x in events)==int(expect["candidate_calls"])
            assert sum(int(x[5]) for x in events)==int(expect["candidate_exhausted"])
            assert sum(int(x[4]) for x in coords)==int(expect["control_calls"])
            assert sum(int(x[5]) for x in coords)==int(expect["control_exhausted"])
            assert int(np.sum(values<=2))==int(expect["candidate_within2"])
        count=next(csv.DictReader(audit.open()))
        result.append(dict(case=case,driver=driver,solves=len(events),calls=sum(int(x[4]) for x in events),
                      exhausted=sum(int(x[5]) for x in events),raw_sha256=sha(data),**count,**stats))
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    subs=p.add_subparsers(dest="mode",required=True)
    prep=subs.add_parser("prepare");prep.add_argument("source",type=Path)
    run=subs.add_parser("check")
    for name in ("history","work","binary","result"):
        run.add_argument(name,type=Path)
    a=p.parse_args()
    if a.mode=="prepare":
        prepare(a.source.resolve())
    else:
        a.work.mkdir()
        rows=evaluate(a.history.resolve(),a.work.resolve(),a.binary.resolve())
        a.result.write_text(json.dumps(dict(source=CANDIDATE,observer_sha256=sha(a.binary.read_bytes()),libraries=rows),indent=2)+"\n")
