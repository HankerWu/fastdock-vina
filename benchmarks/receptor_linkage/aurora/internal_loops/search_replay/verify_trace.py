"""Read the archived runs independently; do not import the replay/evaluation code."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import tarfile


HERE = Path(__file__).resolve().parent


def unpack_events(data):
    events = defaultdict(dict)
    coordinates, evals = {}, defaultdict(list)
    for line in data.decode().splitlines():
        trajectory, step, stage, payload = line.split(maxsplit=3)
        key = int(trajectory), int(step)
        if stage == "path":
            kind, payload = payload.split()
            stage = kind+"_path"
        if stage == "eval":
            kind, payload = payload.split(maxsplit=1)
            evals[(*key,kind)].append(payload)
        elif stage == "end_xyz":
            kind, payload = payload.split(maxsplit=1)
            coordinates[(*key,kind)] = list(map(float,payload.split()))
        else:
            assert stage not in events[key]
            events[key][stage] = payload
    return events, coordinates, evals


def verify(work):
    digest = lambda b: hashlib.sha256(b).hexdigest()
    runs = list(csv.DictReader((HERE/"results/runs.csv").open()))
    findings = list(csv.DictReader((HERE/"results/divergences.csv").open()))
    order = ("start","proposal","hunt_path","hunt","accept","second_path","second","saved")
    costs = defaultdict(int)
    summary = defaultdict(int)
    with tarfile.open(HERE/"traces.tar.xz") as tf:
        # Read each member once; compressed tar seeks otherwise decompress repeatedly.
        raw = {Path(m.name).name: tf.extractfile(m).read() for m in tf if m.isfile()}
    for r in runs:
        drug, seed, arm, phase = (r[k] for k in ("drug","seed","arm","phase"))
        tag = "_".join((drug,seed,arm,phase))
        trace, pose, log = (raw[tag+ext] for ext in (".trace",".pdbqt",".log"))
        assert digest(trace)==r["trace_sha256"] and digest(pose)==r["pose_sha256"]
        path, ref = ("templates","1OL7_released") if arm == "released" else ("base","dry")
        archived = work/path/"blind"/f"{ref}_{drug}_{seed}.pdbqt"
        assert pose == archived.read_bytes()
        select = lambda b: [l for l in b.splitlines() if l.startswith(
            (b"Search trajectory:",b"Hunt chain:",b"Pocket start:",b"Search counters:"))]
        assert select(log)==select(archived.with_suffix(".log").read_bytes()) and len(select(log))==7
        costs[phase] += sum(map(int,re.findall(rb"\b(?:hunt|second|explicit)_evals=(\d+)",log)))
        if arm == "far":
            assert trace == raw[f"{drug}_{seed}_dry_stages.trace"]
            summary["far_equals_dry"] += 1
        summary["runs_match_archived_producer"] += 1
    for r in findings:
        tr, prefix = int(r["trajectory"]), f"{r['drug']}_{r['seed']}"
        (a, ax, ae), (b, bx, be) = [unpack_events(raw[prefix+f"_{arm}_focus.trace"])
                                      for arm in ("dry","released")]
        changed = [(s, i, stage) for t,s in a.keys()|b.keys() if t==tr
                   for i,stage in enumerate(order) if a[t,s].get(stage)!=b[t,s].get(stage)]
        step, _, stage = min(changed)
        assert step==int(r["first_step"]) and stage==r["first_stage"]
        kind = stage.removesuffix("_path")
        key = tr,step,kind
        for i,(left,right) in enumerate(zip(ae[key],be[key]),1):
            if left == right:
                continue
            l_head, l_xyz = left.split(" xyz ")
            r_head, r_xyz = right.split(" xyz ")
            l_head, l_grad = l_head.split(" grad ")
            r_head, r_grad = r_head.split(" grad ")
            l_no, l_energy, l_conf = l_head.split(maxsplit=2)
            r_no, r_energy, r_conf = r_head.split(maxsplit=2)
            assert i==int(l_no)==int(r_no)==int(r["first_changed_eval"])
            assert l_conf==r_conf and l_xyz==r_xyz and l_grad!=r_grad
            assert float(r_energy)-float(l_energy)==float(r["eval_energy_delta"])
            break
        else:
            raise ValueError(f"Failed to find a changed evaluation: {key}")
        d = math.dist(ax[key], bx[key]) / math.sqrt(len(ax[key])/3)
        assert math.isclose(d,float(r["local_end_same_atom_rms_A"]),abs_tol=1e-14)
        summary[stage] += 1
        summary["identical_first_local_return_geometry"] += int(d==0)
        different = [(s, j, stage) for t,s in a.keys() & b.keys() if t==tr
                     for j,stage in enumerate(("hunt","second"))
                     if stage in a[t,s] and stage in b[t,s]
                     and a[t,s][stage].split(maxsplit=2)[2] != b[t,s][stage].split(maxsplit=2)[2]]
        if different:
            s, _, stage = min(different)
            assert s==int(r["first_changed_local_geometry_step"])
            key = tr,s,stage
            d = math.dist(ax[key],bx[key])/math.sqrt(len(ax[key])/3)
            assert math.isclose(d,float(r["first_changed_local_geometry_rms_A"]),abs_tol=1e-14)
            assert s<int(r["accept"])
        else:
            assert not r["first_changed_local_geometry_step"] and not r["accept"]
            summary["no_changed_return_geometry"] += 1
    return {"counts":dict(summary),"search_and_explicit_evaluations":dict(costs)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    print(json.dumps(verify(parser.parse_args().work),indent=2,sort_keys=True))
