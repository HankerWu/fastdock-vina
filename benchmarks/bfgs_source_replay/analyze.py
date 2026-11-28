"""Analyze native and counterfactual records with the frozen receptor-frame symmetry maps."""
import argparse
from collections import defaultdict
from itertools import groupby
import json
from pathlib import Path
import statistics
import numpy as np

from run_replay import ARMS, CASES, run_panel


def native(path):
    trajectories, explicit = [], defaultdict(list)
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        typ = fields[0]
        if typ == "T":
            trajectories.append([])
        elif typ in {"I", "O", "F"}:
            explicit[typ].append(fields)
        else:
            if typ == "S":
                trajectories[-1].append(defaultdict(list))
            step = trajectories[-1][-1]
            assert int(fields[1]) == len(trajectories[-1]) - 1
            step[typ].append(fields)
    return trajectories, explicit


def best(rows):
    return min(rows, default=float("nan"))


def main(work, output):
    output.mkdir(parents=True, exist_ok=True)
    tables = defaultdict(list)
    costs = json.loads((work / "costs.json").read_text())
    for query, receiver, seed in sorted(CASES):
        case = f"{query}-{receiver}-{seed}"
        path = work / case
        fx = work / "fixtures/reciprocal"
        q, r = [json.loads((fx / p / "preparation.json").read_text()) for p in (query, receiver)]
        ref = run_panel.replay.reciprocal.reference(q, r)
        # Trace coordinate order is ligand.cont PDBQT input order, excluding H.
        serials = [int(l[6:11]) for l in (fx / query / "ligand.pdbqt").read_text().splitlines()
                   if l.startswith(("ATOM  ", "HETATM")) and l[77:].strip() not in {"H", "HD", "HS"}]
        assert set(serials) == set(q["serial"]) and len(serials) == len(q["serial"])
        order = [serials.index(i) for i in q["serial"]]

        def rmsd(fields):
            xyz = np.fromstring(fields[8], dtype=int, sep=" ").reshape(-1, 3)[order] / 1000
            return float(np.mean(np.sum((xyz[None, :, :] - ref)**2, axis=-1), axis=1).min()**0.5)

        runs = {arm: native(path / (arm + ".native")) for arm in ARMS}
        for tr, (old, new) in enumerate(zip(runs[ARMS[0]][0], runs[ARMS[1]][0]), 1):
            def first(which):
                for step, (a, b) in enumerate(zip(old, new)):
                    for phase in ("S", "P", "H", "A", "D", "M"):
                        if which == "cost":
                            if phase not in ("D", "H") or not a[phase] or not b[phase]:
                                continue
                            x, y = a[phase][0][4], b[phase][0][4]
                        else:
                            # P carries the reused output buffer's stale energy; S at zero
                            # starts from a randomized pose with zero placeholder energy.
                            cols = slice(3, 4) if phase == "P" else slice(2, 4)
                            x = [v[cols] for v in a[phase]]
                            y = [v[cols] for v in b[phase]]
                        if x != y:
                            return step, phase
                return -1, "none"

            split, stage = first("state")
            csplit, cstage = first("cost")
            row = dict(case=case, trajectory=tr, state_step=split, state_stage=stage,
                       cost_step=csplit, cost_stage=cstage)
            if split >= 0:
                o, n = old[split], new[split]
                row.update(same_source=o["S"][0][2:4] == n["S"][0][2:4],
                           same_proposal=o["P"][0][3] == n["P"][0][3])
                for idx, event in enumerate((o, n)):
                    h = event["H"][0]
                    k = ("old", "repaired")[idx]
                    row.update({k+"_hunt_evals":int(h[4]), k+"_hunt_e": float(h[2]),
                                k+"_hunt_rmsd":rmsd(h), k+"_accepted": bool(event["A"]),
                                k+"_saved": bool(event["D"]), k+"_first_exhausted_angle":float(h[7]),
                                k+"_non_descent":int(h[6])})
                    stage_end = event[stage][0]
                    row.update({k+"_split_e":float(stage_end[2]), k+"_split_evals":int(stage_end[4]),
                                k+"_split_rmsd":rmsd(stage_end)})
            tables["first_divergence"].append(row)

        for arm in ARMS:
            trajectories, explicit = runs[arm]
            cost = next(s for s in costs if s["case"]==case and s["kind"]=="native" and s["driver"]==arm)
            released = run_panel.replay.reciprocal.measure(path/(cost["key"]+".pdbqt"), q["serial"], ref)
            # The human-facing mode table has no energy-window cutoff. The actual
            # file writer additionally applies the frozen 100 kcal/mol window.
            displayed = explicit["F"]
            published = [e for e in displayed if float(e[2]) <= float(displayed[0][2]) + 100]
            assert len(released) == len(published)
            assert max(abs(rmsd(e)-val) for e,val in zip(published, released)) < 0.001
            raw_models = (path/(cost["key"]+".pdbqt")).read_text().split("ENDMDL")[:-1]
            for e, model in zip(published, raw_models):
                atom_xyz = {int(l[6:11]):[float(l[k:k+8]) for k in (30,38,46)]
                            for l in model.splitlines() if l.startswith(("ATOM  ","HETATM"))}
                xyz = np.fromstring(e[8], dtype=int, sep=" ").reshape(-1,3) / 1000
                assert np.max(np.abs(xyz-np.asarray([atom_xyz[n] for n in serials]))) < 0.00101
            inputs = {e[3]:int(e[1]) for e in explicit["I"]}
            outputs = {e[3]:int(e[1]) for e in explicit["O"]}
            assert len(inputs) == len(outputs) == len(explicit["I"]) == len(explicit["O"])
            origins = {}
            flat = [step for traj in trajectories for step in traj]
            for tr, traj in enumerate(trajectories, 1):
                mate = next(k for k in ARMS if k != arm)
                length = min(len(traj), len(runs[mate][0][tr-1]))
                seq = [rmsd(s["H"][0]) for s in traj]
                accepted = [x for x,s in zip(seq,traj) if s["A"]]
                saved, retained_ever, per_step = [], [], {}
                previous, final_score_errors = set(), {}
                n_insertions, n_evictions = 0, 0
                for i, event in enumerate(traj):
                    if not event["D"]:
                        continue
                    val = event["D"][0]
                    state = val[3]
                    value = rmsd(val)
                    saved.append(value)
                    current = {e[3] for e in event["M"]}
                    assert len(current) == len(event["M"]) and current - previous <= {state}
                    if state in current - previous:
                        retained_ever.append(value)
                        final_score_errors[state] = value
                        per_step[state] = i
                        n_insertions += 1
                    n_evictions += len(previous - current)
                    previous = current
                for st in previous & inputs.keys():
                    assert st not in origins
                    origins[st] = (tr, per_step[st])
                final = [final_score_errors[s] for s in previous]
                runs_near = [sum(1 for _ in g) for near, g in groupby(x <= 2 for x in seq) if near]
                row = dict(case=case, driver=arm, trajectory=tr, hunts=len(seq), prefix_n=length,
                           prefix_best=best(seq[:length]), prefix_near=sum(x<=2 for x in seq[:length]),
                           hunt_best=best(seq), hunt_near=sum(x<=2 for x in seq),
                           longest_near_run=max(runs_near,default=0), accepted=len(accepted),
                           accepted_best=best(accepted), accepted_near=sum(x<=2 for x in accepted),
                           saved_attempts=len(saved), saved_best=best(saved), saved_near=sum(x<=2 for x in saved),
                           inserted=n_insertions, evicted=n_evictions, inserted_best=best(retained_ever),
                           inserted_near=sum(x<=2 for x in retained_ever), final_pool=len(final),
                           final_best=best(final), final_near=sum(x<=2 for x in final))
                tables["path_attrition"].append(row)
            assert origins.keys() == inputs.keys()
            by_index = {i:(state, *origins[state]) for state,i in inputs.items()}
            ends = {int(e[1]):e for e in explicit["O"]}
            for e in published:
                idx=outputs[e[3]]
                _,traj,step = by_index[idx]
                saved = trajectories[traj-1][step]["D"][0]
                tables["published_lineage"].append(dict(case=case,driver=arm,rank=int(e[1])+1,
                    trajectory=traj,step=step,merged_index=idx,map_energy=float(saved[2]),
                    saved_rmsd=rmsd(saved),refined_rmsd=rmsd(ends[idx]),published_rmsd=released[int(e[1])],
                    score=float(e[2])))
            by_solver = {}
            for solver in ARMS:
                raw = (path/f"{arm}-{solver}.shadow").read_text().splitlines()
                rows = [l.split("\t") for l in raw]
                assert len(rows)==len(flat) and all(int(r[1])==i for i,r in enumerate(rows))
                if solver == arm:
                    assert all(r[2:]==s["H"][0][2:] for r,s in zip(rows,flat))
                by_solver[solver]=rows
            old, repaired = [by_solver[a] for a in ARMS]
            denergy=[float(b[2])-float(a[2]) for a,b in zip(old,repaired)]
            dr=[rmsd(b)-rmsd(a) for a,b in zip(old,repaired)]
            oldr,newr=[list(map(rmsd,rows)) for rows in (old,repaired)]
            tables["locked_proposals"].append(dict(case=case,driver=arm,n=len(flat),
                old_best=best(oldr),repaired_best=best(newr),
                old_near=sum(x<=2 for x in oldr),repaired_near=sum(x<=2 for x in newr),
                n_equal_state=sum(a[3]==b[3] for a,b in zip(old,repaired)),
                same_solver_reproduced=len(flat),delta_energy_median=statistics.median(denergy),
                delta_rmsd_median=statistics.median(dr),energy_worse=sum(v>1e-7 for v in denergy),
                energy_better=sum(v< -1e-7 for v in denergy),rmsd_worse=sum(v>0.001 for v in dr),
                rmsd_better=sum(v< -0.001 for v in dr),joint_worse=sum(v>1e-7 and w>0.001 for v,w in zip(denergy,dr)),
                old_evals=sum(int(a[4]) for a in old),repaired_evals=sum(int(b[4]) for b in repaired),
                repaired_exhausted=sum(int(b[5])>0 for b in repaired),
                repaired_non_descent=sum(int(b[6])>0 for b in repaired),
                exhausted_last_angle_gt_pi=sum(float(b[7])/512>np.pi for b in repaired)))
    for name,rows in tables.items():
        run_panel.csvfile(output/(name+".csv"),rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("output",type=Path)
    args=p.parse_args()
    main(args.work,args.output)
