"""Independently audit released coordinates, locked returns and phase attrition."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"bfgs_source_replay"))
from verify import reference

def check(work,tables):
    def rows(table):
        return list(csv.DictReader((tables/(table+".csv")).open()))
    def records(path):
        return [line.split("\t") for line in path.read_text().splitlines()]
    checked={"published":0,"paths":0,"same_solver":0,"locked_returns":0,"corners":0}
    largest=0
    pub=rows("published")
    locked=rows("locked_solves")
    attrition=rows("attrition")
    archive=list(csv.DictReader((Path(__file__).resolve().parents[1]/"search_trial_panel/results/cells.csv").open()))
    for cell in rows("free_cells"):
        case,mode=cell["case"],cell["mode"]
        q,r,seed=case.split("-")
        inputs=work/"fixtures/reciprocal"
        query,receiver=(json.loads((inputs/p/"preparation.json").read_text()) for p in (q,r))
        mapped=reference(query,receiver)
        xyz_order=[int(l[6:11]) for l in (inputs/q/"ligand.pdbqt").read_text().splitlines()
                if l.startswith(("ATOM  ","HETATM")) and int(l[6:11]) in query["serial"]]
        reorder=[xyz_order.index(n) for n in query["serial"]]
        def measure(coords):
            diff=coords-mapped
            return np.sqrt(np.einsum("sij,sij->s",diff,diff).min()/len(reorder))
        def event_rmsd(rec):
            return measure(np.fromstring(rec[8],sep=" ").reshape(-1,3)[reorder]*0.001)
        def compare(x,y):
            nonlocal largest
            d=abs(x-float(y))
            assert d<1.e-9
            largest=max(largest,d)
        dest=work/case/mode
        pool,=dest.glob("*.pdbqt")
        if mode in ("00","11"):
            a=next(s for s in archive if s["key"]==pool.stem)
            assert hashlib.sha256(pool.read_bytes()).hexdigest()==a["pool_sha256"]
            checked["corners"]+=1
        for rank,block in enumerate(pool.read_text().split("ENDMDL")[:-1],1):
            coords={int(l[6:11]): [float(l[k:k+8]) for k in (30,38,46)]
                    for l in block.splitlines() if l.startswith(("ATOM  ","HETATM"))}
            actual=measure(np.array([coords[n] for n in query["serial"]]))
            expected=next(x for x in pub if x["mode"]==mode and x["case"]==case and int(x["rank"])==rank)
            compare(actual,expected["top_rmsd"])
            checked["published"]+=1
        trajectories=[]
        free=records(dest/"free.native")
        for f in free:
            if f[0]=="T": trajectories.append([])
            elif f[0] in ("A","H","D","M"):trajectories[-1].append(f)
        for i,events in enumerate(trajectories,1):
            a=next(x for x in attrition if x["case"]==case and x["mode"]==mode and int(x["trajectory"])==i)
            for event,name in (("H","hunt"),("D","saved")):
                xs=[event_rmsd(f) for f in events if f[0]==event]
                assert len(xs)==int(a[name]) and sum(x<=2 for x in xs)==int(a[name+"_near"])
                compare(min(xs),a[name+"_best"])
            assert sum(f[0]=="A" for f in events)==int(a["accepted"])
            checked["paths"]+=1
        if mode in ("00","11"):
            raw=records(dest/"locked.trace")
            if mode=="00":
                expected=[f for f in free if f[0] in ("H","D")]
                assert len(raw)==len(expected) and all(x[0]==y[0] and x[2:-1]==y[2:-1] for x,y in zip(raw,expected))
                checked["same_solver"]+=len(raw)
            for phase in ("H","D"):
                values=[x for x in raw if x[0]==phase]
                row=next(s for s in locked if s["case"]==case and s["mode"]==mode and s["phase"]==phase)
                compare(min(event_rmsd(x) for x in values),row["best"])
                assert sum(int(x[4]) for x in values)==int(row["evals"])
                assert sum(int(x[13]) for x in values)==int(row["directions_with_clipped_translation"])
                checked["locked_returns"]+=len(values)
    return dict(**checked,max_abs_geometry_discrepancy_A=largest)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("tables",type=Path)
    args=p.parse_args()
    print(json.dumps(check(args.work,args.tables),indent=2))
