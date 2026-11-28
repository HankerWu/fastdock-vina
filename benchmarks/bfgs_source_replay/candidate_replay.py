"""Run and compare the two MC-configured grid optimizers on every frozen proposal."""
import argparse
import json
import os
from pathlib import Path
import numpy as np

from run_replay import ARMS, CASES, digest, run_panel


def compare(work, output, binaries=None):
    results, costs = [], []
    fixtures = work / "fixtures/reciprocal"
    for query, receiver, seed in sorted(CASES):
        case = f"{query}-{receiver}-{seed}"
        path = work / case
        q,r = [json.loads((fixtures / name / "preparation.json").read_text()) for name in (query,receiver)]
        reference = run_panel.replay.reciprocal.reference(q,r)
        serials = [int(l[6:11]) for l in (fixtures/query/"ligand.pdbqt").read_text().splitlines()
                   if l.startswith(("ATOM  ", "HETATM")) and l[77:].strip() not in {"H","HD","HS"}]
        order = [serials.index(i) for i in q["serial"]]
        for driver in ARMS:
            returns={}
            for engine in ("ed0d42d","0bb827e"):
                prefix=f"{driver}-{engine}"
                shadow=path/(prefix+".shadow")
                if binaries:
                    binary=binaries[engine]
                    os.environ["FD_BFGS_SHADOW"]=str(path/(driver+".inputs"))
                    os.environ["FD_BFGS_TRACE"]=str(shadow)
                    cmd=[str(binary),"--ligand",str(fixtures/query/"ligand.pdbqt"),"--receptor",
                         str(fixtures/receiver/"receptor.pdbqt"),"--cpu","1","--seed",str(seed),
                         "--max_evals","12000","--exhaustiveness","3","--num_modes","9",
                         "--energy_range","100","--spacing","0.375","--out",str(path/(prefix+".unused.pdbqt"))]
                    for axis,value in zip("xyz",r["center"]):
                        cmd+=["--size_"+axis,"24","--center_"+axis,str(value)]
                    stats=run_panel.execute(cmd,path/(prefix+".log"),180)
                    assert stats["status"]=="ok"
                    costs.append(dict(case=case,driver=driver,evaluator=engine,kind="candidate_shadow",
                                **stats,trace_sha256=digest(shadow),binary_sha256=digest(binary)))
                raw=[l.split("\t") for l in shadow.read_text().splitlines()]
                assert raw and all(int(f[1])==i for i,f in enumerate(raw))
                def error(row):
                    xyz=np.fromstring(row[8],dtype=int,sep=" ").reshape(-1,3)[order]/1000
                    return float(np.sqrt(((xyz[None]-reference)**2).sum(axis=-1).mean(axis=1)).min())
                returns[engine]=(raw,list(map(error,raw)))
            (a,ra),(b,rb)=returns.values()
            assert len(a)==len(b)
            # f0 and directional slope -||g||^2 for unit initial inverse Hessian.
            assert all(x[9:11]==y[9:11] for x,y in zip(a,b))
            de=np.array([float(y[2])-float(x[2]) for x,y in zip(a,b)])
            dr=np.array(rb)-np.array(ra)
            row=dict(case=case,driver=driver,n=len(a),equal_initial_E_and_pg=len(a),
                     median_paired_energy=float(np.median(de)),median_paired_rmsd=float(np.median(dr)),
                     lower_energy=int((de < -1e-7).sum()),higher_energy=int((de>1e-7).sum()),
                     better_rmsd=int((dr < -0.001).sum()),worse_rmsd=int((dr>0.001).sum()))
            for key, (vals, rr) in zip(("stable","candidate"),returns.values()):
                exhausted=[v for v in vals if int(v[5])]
                stops=[v for v in vals if int(v[6])]
                row.update({key+"_exhausted":len(exhausted),key+"_non_descent":len(stops),
                    key+"_non_descent_positive":sum(float(v[11])>0 for v in stops),
                    key+"_exhausted_nonnegative_pg":sum(float(v[11])>=0 for v in exhausted),
                    key+"_last_angle_gt_pi":sum(float(v[12])>np.pi for v in exhausted),
                    key+"_max_last_angle":max((float(v[12]) for v in exhausted),default=0),
                    key+"_evals":sum(int(v[4]) for v in vals),key+"_best_rmsd":min(rr),
                    key+"_near":sum(x<=2 for x in rr)})
            results.append(row)
    run_panel.csvfile(output/"candidate_locked.csv",results)
    if binaries:
        (work/"candidate_costs.json").write_text(json.dumps(costs,indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("--binaries",nargs=2,type=Path,metavar=("STABLE","CANDIDATE"))
    a=p.parse_args()
    compare(a.work.resolve(),a.output,dict(zip(("ed0d42d","0bb827e"),map(Path.resolve,a.binaries))) if a.binaries else None)
