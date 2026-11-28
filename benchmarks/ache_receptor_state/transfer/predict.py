"""Transfer every blind donor pose, optimize and publish using recipient scores alone."""
import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "blind"))
from dock import RECEIVERS, SEEDS, load_inputs, save_json, sha, table, tick
sys.path.insert(0, str(HERE.parent / "diagnostic"))
from measure import move, model_positions, score

BLIND = HERE.parent / "blind/results"


def rms(x, y):
    return float(np.sqrt(np.mean(np.sum((x-y)**2, axis=1))))


def run(binary, args, job, name):
    start = tick()
    with (job / (name + ".log")).open("wb") as stream:
        cmd = subprocess.run([str(binary), *args], cwd=job, stdout=stream, stderr=subprocess.STDOUT)
    wall, cpu = tick()-start
    return dict(task=name, exit_code=cmd.returncode, wall_s=wall, cpu_s=cpu)


def inspect_one(rec, source, seed, mode, text, template, matrix, delta, center, heavy, binary, observer, work):
    job = work / "raw" / rec / f"{source}_{seed}_{mode:02d}"
    job.mkdir()
    clock = tick()
    sites = model_positions(text.splitlines())
    moved = move(template, sites, matrix, delta)
    (job / "start.pdbqt").write_text(moved)
    starts = model_positions(moved.splitlines())
    if set(starts) != set(sites):
        raise ValueError("Missing ligand atom in transformed start")
    frame = ["--receptor", "../../"+rec+".pdbqt", "--cpu", "1", "--seed", str(seed)]
    for axis, value in zip("xyz", center, strict=True):
        frame += ["--center_"+axis, str(value), "--size_"+axis, "36"]
    local = frame+["--ligand","start.pdbqt","--local_only","--out","end.pdbqt"]
    final = frame+["--ligand","end.pdbqt","--score_only","--unbound_energy","0"]
    observed = local[:-1]+["observer.pdbqt"]
    save_json(job/"commands.json", dict(local=local, score=final, observer=observed))
    record = dict(receiver=rec, source=source, seed=seed, mode=mode, status="ok",
                  start_sha256=sha(moved.encode()), end_sha256=None, observer_match=None,
                  box_max_initial_A=float(np.abs(np.stack(list(starts.values()))-center).max()),
                  box_max_final_A=None, score=None, inter=None, intra=None, torsion=None, unbound=None,
                  objective_in=None, objective_out=None, gradient_calls=None, local_shift_A=None,
                  production_wall_s=None, production_cpu_s=None, observer_wall_s=None, observer_cpu_s=None)
    timings, coords = [], None
    for name, args in (("local", local), ("score", final)):
        timing = run(binary, args, job, name)
        timings.append(timing)
        if timing["exit_code"]:
            record["status"] = name + "_failed"
            break
    if record["status"] == "ok":
        ends = model_positions((job / "end.pdbqt").read_text().splitlines())
        if ends.keys() != starts.keys():
            raise ValueError("Missing ligand atom after native optimization")
        coords = np.stack([ends[i] for i in heavy])
        record.update(dict(zip(("score","inter","intra","torsion","unbound"), score(job/"score.log"), strict=True)))
        record.update(end_sha256=sha((job/"end.pdbqt").read_bytes()),
                      box_max_final_A=float(np.abs(np.stack(list(ends.values()))-center).max()),
                      local_shift_A=rms(coords, np.stack([starts[i] for i in heavy])))
    record["production_wall_s"], record["production_cpu_s"] = tick()-clock
    if record["status"] == "ok":
        timing = run(observer, observed, job, "observer")
        timings.append(timing)
        observed_file = job/"observer.pdbqt"
        record["observer_wall_s"], record["observer_cpu_s"] = timing["wall_s"], timing["cpu_s"]
        record["observer_match"] = observed_file.is_file() and sha(observed_file.read_bytes()) == record["end_sha256"]
        obs = re.findall(r"REFINE_OBSERVE\s+(\S+)\s+(\S+)\s+(\d+)", (job/"observer.log").read_text())
        if record["observer_match"] and obs and timing["exit_code"] == 0:
            record.update(objective_in=float(obs[0][0]), objective_out=float(obs[-1][1]),
                          gradient_calls=sum(int(x[2]) for x in obs))
    save_json(job/"timing.json", timings)
    return record, coords


def pool(blind, rec, arm, seed, runs):
    path = blind / "raw" / rec / f"{arm}_{seed}" / "dry.pdbqt"
    data = path.read_bytes()
    if sha(data) != runs[rec, arm, seed]["output_sha256"]:
        raise ValueError("Source dry pool digest changed")
    blocks = data.decode().split("ENDMDL")
    if blocks[-1].strip():
        raise ValueError("Unterminated published PDBQT")
    return blocks[:-1]


def select(origins, jobs, coords):
    valid = [j for j in origins if jobs[j]["status"] == "ok"]
    ordered = sorted(valid, key=lambda j: jobs[j]["score"])
    accepted, results = [], []
    for j in ordered:
        score_value = jobs[j]["score"]
        duplicate = None
        for k in accepted:
            if rms(coords[j], coords[k]) < 1:
                duplicate = k
                break
        if score_value > jobs[ordered[0]]["score"] + 100:
            status = "window"
        elif duplicate is not None:
            status = "duplicated"
        elif len(accepted) >= 9:
            status = "rank_cut"
        else:
            accepted.append(j)
            status = "retained"
        results.append(dict(source=j[1], source_seed=j[2], source_mode=j[3], score=score_value,
                            decision=status, rank=len(accepted) if status=="retained" else None,
                            duplicate_of=str(duplicate) if duplicate is not None else None))
    for j in origins:
        if jobs[j]["status"] != "ok":
            results.append(dict(source=j[1], source_seed=j[2], source_mode=j[3], score=None,
                                decision="failed", rank=None, duplicate_of=None))
    return results


def main(work, out, binary, observer, blind):
    work.mkdir(parents=True, exist_ok=False)
    out.mkdir(parents=True, exist_ok=False)
    start = tick()
    prep = load_inputs(work)
    template = (work/"HUP.pdbqt").read_text()
    heavy = [int(l[6:11]) for l in template.splitlines() if l.startswith(("ATOM","HETATM")) and l[77:].strip() not in ("H","HD")]
    if len(heavy) != 18:
        raise ValueError("Expected 18 HUP heavy atoms")
    runs = {(r["receiver"], r["arm"], int(r["seed"])):r for r in csv.DictReader((blind/"runs.csv").open())}
    r0, r1 = RECEIVERS
    rot = np.array(prep["noncognate_fit"]["rotation"])
    trans = np.array(prep["noncognate_fit"]["translation_A"])
    source_hashes = {}
    for rec in RECEIVERS:
        path = work/"raw"/rec
        path.mkdir(parents=True)
        shutil.copyfile(work/f"{rec}.pdbqt",work/"raw"/f"{rec}.pdbqt")
    jobs, coords, numbers = {}, {}, {}
    for rec in RECEIVERS:
        center = np.array(prep["centers_A"][rec])
        for source in RECEIVERS:
            if rec == source:
                matrix, delta = np.eye(3), np.zeros(3)
            elif source == r0:
                matrix, delta = rot, trans
            else:
                matrix, delta = rot.T, -trans @ rot.T
            for seed in SEEDS:
                blocks = pool(blind, source, "pro24k", seed, runs)
                source_hashes[f"{source}/pro24k_{seed}"] = runs[source,"pro24k",seed]["output_sha256"]
                numbers[source,seed] = len(blocks)
                for mode, block in enumerate(blocks, 1):
                    key = rec, source, seed, mode
                    jobs[key], coords[key] = inspect_one(*key, block, template, matrix, delta, center, heavy,
                                                         binary, observer, work)
    table(out/"jobs.csv", list(jobs.values()))
    rows, costs = [], []
    for rec in RECEIVERS:
        other = r1 if rec==r0 else r0
        for i, seed in enumerate(SEEDS):
            next_seed = SEEDS[(i+1) % len(SEEDS)]
            for arm, pools in (("recipient",[(rec,seed)]), ("cross",[(rec,seed),(other,seed)]),
                               ("same",[(rec,seed),(rec,next_seed)]), ("double",[(rec,seed)])):
                arm_clock=tick()
                source_arm = "pro48k" if arm=="double" else "pro24k"
                used_runs = [runs[s,source_arm,t] for s,t in pools]
                keys = [(rec,s,t,m) for s,t in pools for m in range(1,numbers[s,t]+1)] if arm!="double" else []
                if arm=="double":
                    blocks = pool(blind, rec, source_arm, seed, runs)
                    source_hashes[f"{rec}/{source_arm}_{seed}"]=used_runs[0]["output_sha256"]
                    items = [dict(source=rec, source_seed=seed, source_mode=m, score=float(re.search(
                        r"REMARK VINA RESULT:\s*([-\d.]+)",b)[1]), decision="retained",rank=m,duplicate_of=None)
                             for m,b in enumerate(blocks,1)]
                else:
                    items = select(keys, jobs, coords)
                select_s = tick()-arm_clock
                rows.extend(dict(receiver=rec,seed=seed,arm=arm,**item) for item in items)
                search_time=np.sum([[float(run[k]) for k in ("wall_s","cpu_s")] for run in used_runs],axis=0)
                refinement_time=np.sum([[jobs[j][k] for k in ("production_wall_s","production_cpu_s")] for j in keys],
                                       axis=0) if keys else np.zeros(2)
                observer_time=np.sum([[jobs[j][k] or 0 for k in ("observer_wall_s","observer_cpu_s")] for j in keys],
                                     axis=0) if keys else np.zeros(2)
                total=search_time+refinement_time+select_s
                evals=sum(int(run[k]) for run in used_runs for k in ("hunt_evals","second_evals","explicit_evals"))
                counted=all(jobs[j]["gradient_calls"] is not None for j in keys)
                local_calls=sum(jobs[j]["gradient_calls"] for j in keys) if counted else None
                costs.append(dict(receiver=rec,seed=seed,arm=arm,searches=len(pools),search_wall_s=search_time[0],
                                  search_cpu_s=search_time[1],search_gradient_evals=evals,local_starts=len(keys),
                                  local_gradient_evals=local_calls,final_score_calls=sum(jobs[j]["status"]=="ok" for j in keys),
                                  local_and_score_wall_s=refinement_time[0],local_and_score_cpu_s=refinement_time[1],
                                  select_wall_s=select_s[0],select_cpu_s=select_s[1],
                                  charged_wall_s=total[0],charged_cpu_s=total[1],cells_per_cpu_s=1/total[1],
                                  audit_wall_s=observer_time[0],audit_cpu_s=observer_time[1]))
    table(out/"selections.csv",rows)
    table(out/"costs.csv",costs)
    with tarfile.open(out/"raw.tar.xz","w:xz") as tar:
        for rec in RECEIVERS:
            tar.add(work/"raw"/rec, arcname=rec)
    save_json(out/"receipt.json", dict(binary_sha256=sha(binary.read_bytes()),observer_sha256=sha(observer.read_bytes()),
                                       inputs=prep,pool_sha256=source_hashes,raw_sha256=sha((out/"raw.tar.xz").read_bytes()),
                                       n_jobs=len(jobs),failures=[list(k) for k,r in jobs.items() if r["status"]!="ok"],
                                       total_wall_cpu_s=(tick()-start).tolist()))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("work","out","pro","observer"):
        p.add_argument(arg,type=Path)
    p.add_argument("--blind",type=Path,default=BLIND)
    a=p.parse_args()
    main(a.work.resolve(),a.out.resolve(),a.pro.resolve(),a.observer.resolve(),a.blind.resolve())
