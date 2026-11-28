"""Rebuild intervention effects, phase counts, per-trajectory attrition and pose lineage."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import statistics
import numpy as np
from run_stages import MODES, cases, native, replay_inputs, run_panel, HERE, digest

def best(errors):
    return min(errors,default=None)


def observe_trajectory(traj, error):
    retained = {}
    previous = set()
    inserted, evicted, near_evicted = [], 0, 0
    for i,step in enumerate(traj):
        if not step["D"]:
            continue
        d = step["D"][0]
        current = {m[3] for m in step["M"]}
        new = current - previous
        assert new <= {d[3]}
        if new:
            value = error(d)
            retained[d[3]] = (i,value)
            inserted.append(value)
        evicted += len(previous-current)
        near_evicted += sum(retained[s][1] <= 2 for s in previous-current)
        previous=current
    return {s:retained[s] for s in previous}, inserted, evicted, near_evicted


def divergence(ref, test):
    for i,(a,b) in enumerate(zip(ref,test)):
        for phase in ("S","P","H","A","D","M"):
            cols=slice(3,4) if phase=="P" else slice(2,4)
            if [r[cols] for r in a[phase]] != [r[cols] for r in b[phase]]:
                return i,phase,a["S"][0][2:4]==b["S"][0][2:4],a["P"][0][3]==b["P"][0][3]
    return min(len(ref),len(test)),"length" if len(ref)!=len(test) else "none",None,None


def solve_stats(records):
    calls=len(records)
    n_eval=sum(int(x[4]) for x in records)
    cpu=sum(int(x[-1]) for x in records)/1e9
    return dict(solves=calls,evals=n_eval,evals_per_solve=n_eval/calls,cpu_s=cpu,
          us_per_solve=cpu*1e6/calls,exhausted=sum(int(x[5]) for x in records),
          non_descent=sum(int(x[6]) for x in records),solves_with_clipped_translation=sum(int(r[13])>0 for r in records),
          directions_with_clipped_translation=sum(int(x[13]) for x in records),
          sum_initial_translation_norm_before_A=sum(float(r[14]) for r in records),
          sum_initial_translation_norm_after_A=sum(float(r[15]) for r in records))


def main(work,out):
    out.mkdir(exist_ok=False)
    tables=defaultdict(list)
    costs=json.loads((work/"costs.json").read_text())
    archive=list(csv.DictReader((HERE.parent/"search_trial_panel/results/cells.csv").open()))
    for case, job in cases(work):
        path=work/case
        _,inputs,q,r,ref,seed,cap=job
        serials=[int(l[6:11]) for l in (inputs/q["pdb"]/"ligand.pdbqt").read_text().splitlines()
                 if l.startswith(("ATOM  ","HETATM")) and l[77:].strip() not in {"H","HS","HD"}]
        order=[serials.index(n) for n in q["serial"]]
        def error(rec):
            xyz=np.fromstring(rec[8],dtype=int,sep=" ").reshape(-1,3)[order]/1000
            return float(np.sqrt(np.min(np.mean(np.sum((xyz[None]-ref)**2,axis=-1),axis=1))))

        runs={mode:native(path/mode/"free.native") for mode in MODES}
        for mode,(trajectories,explicit) in runs.items():
            stat=next(row for row in costs if row["case"]==case and row["mode"]==mode and row["phase"]=="free")
            pdbqt=path/mode/(stat["key"]+".pdbqt")
            released=run_panel.replay.reciprocal.measure(pdbqt,q["serial"],ref)
            assert digest(pdbqt)==stat["pool_sha256"]
            if mode in ("00","11"):
                assert digest(pdbqt)==next(row["pool_sha256"] for row in archive if row["key"]==stat["key"])
            displayed=explicit["F"]
            published=[e for e in displayed if float(e[2]) <= float(displayed[0][2])+100]
            assert len(released)==len(published)
            # The F record and PDBQT print their Cartesian coordinates via different rounding routines.
            assert all(abs(error(rec)-value)<0.0011 for rec,value in zip(published,released))
            origins={}
            for t,traj in enumerate(trajectories):
                h=[step["H"][0] for step in traj]
                d=[step["D"][0] for step in traj if step["D"]]
                accept=[error(step["H"][0]) for step in traj if step["A"]]
                final,ins,evicted,near_evicted=observe_trajectory(traj,error)
                for state,(idx,val) in final.items():
                    if state not in origins:
                        origins[state]=(t,idx)
                for label,records in (("hunt",h),("save",d)):
                    tables["trajectory_solves"].append(dict(case=case,mode=mode,trajectory=t+1,phase=label,**solve_stats(records)))
                hr=list(map(error,h)); ds=list(map(error,d))
                final_r=[val for _,val in final.values()]
                prefix=min(len(traj),len(runs["00"][0][t]))
                row=dict(case=case,mode=mode,trajectory=t+1,common_length_with_00=prefix,
                         prefix_hunt_best=best(hr[:prefix]),hunt=len(hr),hunt_best=best(hr),hunt_near=sum(x<=2 for x in hr),
                         accepted=len(accept),accepted_best=best(accept),accepted_near=sum(x<=2 for x in accept),
                         saved=len(ds),saved_best=best(ds),saved_near=sum(x<=2 for x in ds),
                         inserted=len(ins),inserted_best=best(ins),inserted_near=sum(x<=2 for x in ins),
                         evicted=evicted,evicted_near=near_evicted,retained=len(final_r),
                         retained_best=best(final_r),retained_near=sum(x<=2 for x in final_r))
                tables["attrition"].append(row)
            i_by_state={e[3]:int(e[1]) for e in explicit["I"]}
            from_i={idx:origins[state] for state,idx in i_by_state.items()}
            o_by_state={e[3]:int(e[1]) for e in explicit["O"]}
            assert len(from_i)==len(o_by_state)==len(explicit["I"])
            for rank,(rec,value) in enumerate(zip(published,released),1):
                idx=o_by_state[rec[3]]
                tr,step=from_i[idx]
                event=trajectories[tr][step]
                tables["published"].append(dict(case=case,mode=mode,rank=rank,trajectory=tr+1,step=step,
                      merged_index=idx,hunt=error(event["H"][0]),save=error(event["D"][0]),top_rmsd=value,
                      score=float(rec[2])))
            tables["free_cells"].append(dict(case=case,mode=mode,first=released[0],best5=min(released[:5]),
                  oracle=min(released),regret=released[0]-min(released),modes=len(released),
                  hunts=sum(len(t) for t in trajectories),hunt_evals=stat["hunt_evals"],save_evals=stat["second_evals"],
                  explicit_evals=stat["explicit_evals"],cpu_s=stat["cpu_s"],wall_s=stat["wall_s"],rss_kib=stat["peak_rss_kib"],
                  pool_sha256=digest(pdbqt)))
        for lhs,rhs in (("00","10"),("00","01"),("01","11"),("10","11"),("00","11")):
            for t,(a,b) in enumerate(zip(runs[lhs][0],runs[rhs][0]),1):
                i,phase,ss,sp=divergence(a,b)
                tables["divergence"].append(dict(case=case,reference=lhs,test=rhs,trajectory=t,step=i,phase=phase,
                  source_agrees=ss,proposal_agrees=sp))
        src=list(replay_inputs(path/"00/free.native"))
        assert (path/"stable.states").read_text()=="".join(f"{p} {s}\n" for p,s,_ in src)
        controlled={}
        for mode in ("00","11"):
            events=[s.split("\t") for s in (path/mode/"locked.trace").read_text().splitlines()]
            assert len(events)==len(src)
            for k,(obs,(ph,_,real)) in enumerate(zip(events,src)):
                assert int(obs[1])==k and obs[0]==ph
                assert mode!="00" or obs[2:-1]==real[2:-1]
            controlled[mode]=events
            for phase in ("H","D"):
                subset=[r for r in events if r[0]==phase]
                er=list(map(error,subset))
                tables["locked_solves"].append(dict(case=case,mode=mode,phase=phase,best=best(er),
                         near=sum(x<=2 for x in er),**solve_stats(subset)))
        for ph in ("H","D"):
            a,b=([e for e in controlled[mode] if e[0]==ph] for mode in ("00","11"))
            assert all(x[9:11]==y[9:11] for x,y in zip(a,b))
            dr=[error(y)-error(x) for x,y in zip(a,b)]
            de=[float(y[2])-float(x[2]) for x,y in zip(a,b)]
            tables["locked_effect"].append(dict(case=case,phase=ph,n=len(a),identical_state=sum(x[3]==y[3] for x,y in zip(a,b)),
                 energy_delta_median=statistics.median(de),energy_improved=sum(x< -1e-7 for x in de),
                 energy_worsened=sum(x>1e-7 for x in de),rmsd_delta_median=statistics.median(dr),
                 rmsd_improved=sum(x< -0.001 for x in dr),rmsd_worsened=sum(x>0.001 for x in dr)))
    for name,rows in tables.items():
        run_panel.csvfile(out/(name+".csv"),rows)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path);p.add_argument("output",type=Path)
    args=p.parse_args()
    main(args.work,args.output)
