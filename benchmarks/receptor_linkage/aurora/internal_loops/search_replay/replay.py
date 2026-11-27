"""Replay matched two-trajectory searches and isolate the first changed local call."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import time


HERE = Path(__file__).resolve().parent
SEEDS = (60913, 60943, 61211, 61231)
ARMS = ("dry", "released", "far")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse(path):
    events = {}
    for line in path.read_text().splitlines():
        trajectory, step, stage, values = line.split(maxsplit=3)
        if stage == "path":
            call, values = values.split()
            stage = call+"_path"
        key = int(trajectory), int(step), stage
        if stage not in ("eval", "end_xyz"):
            if key in events:
                raise ValueError(f"Repeated trace key {key}")
            events[key] = values
    return events


def localize(a, b, tr):
    stages = ("start", "proposal", "hunt_path", "hunt", "accept", "second_path", "second", "saved")
    keys = sorted({(step, stages.index(stage), stage)
                   for t, step, stage in a.keys() | b.keys() if t == tr})
    first = {}
    for step, _, stage in keys:
        left, right = a.get((tr, step, stage)), b.get((tr, step, stage))
        if left != right and stage not in first:
            first[stage] = step
        if stage in ("hunt", "second") and left and right:
            ln, le, lc = left.split(maxsplit=2)
            rn, re, rc = right.split(maxsplit=2)
            for name, different in ((stage+"_coordinates", lc != rc), (stage+"_energy", le != re),
                                    (stage+"_calls", ln != rn)):
                if different:
                    first.setdefault(name, step)
    if not first:
        return {"first_stage": "none"}
    first_time, _, first_stage = min((first[s], i, s) for i, s in enumerate(stages) if s in first)
    return {"first_stage": first_stage, "first_step": first_time, **first}


def counters(path):
    return [line for line in path.read_text().splitlines()
            if line.startswith(("Search trajectory:", "Hunt chain:", "Pocket start:", "Search counters:"))]


def rows(path, data):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in data for k in r)))
        writer.writeheader()
        writer.writerows(data)


def read_evals(path, tr, step, phase):
    evals = []
    end = None
    for line in path.read_text().splitlines():
        t, s, event, values = line.split(maxsplit=3)
        if int(t) != tr or int(s) != step:
            continue
        if event not in ("eval", "end_xyz"):
            continue
        kind, values = values.split(maxsplit=1)
        if kind != phase:
            continue
        if event == "eval":
            head, xyz = values.split(" xyz ")
            head, grad = head.split(" grad ")
            index, energy, c = head.split(maxsplit=2)
            evals.append((int(index), float(energy), c, tuple(map(float, grad.split())),
                          tuple(map(float, xyz.split()))))
        elif event == "end_xyz" and end is None:
            end = tuple(map(float, values.split()))
    return evals, end


def norm_delta(a, b):
    return sum((x-y)**2 for x, y in zip(a, b, strict=True))**0.5


def summarize_evals(a, b):
    left, xl = a
    right, xr = b
    summary = {"dry_n_eval": len(left), "released_n_eval": len(right),
               "local_end_same_atom_rms_A": norm_delta(xl, xr) / (len(xl)/3)**0.5}
    for l, r in zip(left, right):
        if l != r:
            n, energy, c, grad, xyz = l
            summary.update(first_changed_eval=n, same_conf=(c == r[2]), same_xyz=(xyz == r[4]),
                           eval_energy_delta=r[1]-energy,
                           max_abs_gradient_delta=max(abs(x-y) for x, y in zip(grad, r[3], strict=True)))
            break
    return summary


def run(work, far, binary, out):
    out.mkdir(parents=True, exist_ok=True)
    raw = out / "raw"
    raw.mkdir(exist_ok=True)
    paths = dict(dry=work/"prepared/2X6E.pdbqt", released=work/"template1OL7/trans_model/released/receptor.pdbqt",
                 far=far)
    center = json.loads((work/"prepared/2X6E.json").read_text())["center"]
    cells, cases = [], []
    for drug in ("at9283", "danusertib"):
        for seed in SEEDS:
            ligand = work/"prepared"/f"{drug}.pdbqt"
            events, hashes, focus = {}, {}, None
            for phase, arms in (("stages", ARMS), ("focus", ("dry", "released"))):
                for arm in arms:
                    tag = f"{drug}_{seed}_{arm}_{phase}"
                    trace, pose, log = [raw/(tag+ext) for ext in (".trace", ".pdbqt", ".log")]
                    env = {**os.environ, "FD_REPLAY_TRACE": str(trace)}
                    env.pop("FD_REPLAY_HUNTS", None)
                    if phase == "focus":
                        env["FD_REPLAY_HUNTS"] = ",".join(str(x) for steps in focus for x in steps)
                    cmd = [binary, "--receptor", paths[arm], "--ligand", ligand]
                    for d, x in zip("xyz", center, strict=True):
                        cmd.extend((f"--center_{d}", str(x), f"--size_{d}", "26"))
                    cmd.extend(("--seed", str(seed), "--exhaustiveness", "2", "--cpu", "1", "--num_modes", "9",
                                "--energy_range", "100", "--max_evals", "25000", "--pose_selection", "late",
                                "--search_counters", "--out", pose))
                    start = time.perf_counter()
                    before = resource.getrusage(resource.RUSAGE_CHILDREN)
                    with log.open("w") as f:
                        subprocess.run(cmd, env=env, stdout=f, stderr=subprocess.STDOUT, check=True)
                    duration = time.perf_counter()-start
                    after = resource.getrusage(resource.RUSAGE_CHILDREN)
                    prefix, prev = ("1OL7_released", "templates") if arm == "released" else ("dry", "base")
                    original = work/prev/"blind"/f"{prefix}_{drug}_{seed}.pdbqt"
                    case = dict(drug=drug, seed=seed, arm=arm, phase=phase, trace_sha256=sha(trace),
                                pose_sha256=sha(pose), original_pool_equal=pose.read_bytes()==original.read_bytes(),
                                original_counters_equal=counters(log)==counters(original.with_suffix(".log")),
                                wall_s=duration,
                                cpu_s=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime)
                    cells.append(case)
                    if not case["original_pool_equal"] or not case["original_counters_equal"]:
                        raise RuntimeError(f"Archived producer mismatch, see {log}")
                    if phase == "stages":
                        events[arm] = parse(trace)
                        hashes[arm] = case["trace_sha256"]
                    elif parse(trace) != events[arm]:
                        raise ValueError(f"Instrumentation changed trajectory in {trace}")
                if phase == "stages":
                    if hashes["far"] != hashes["dry"]:
                        raise ValueError(f"Far trace diverged for {drug} seed {seed}")
                    pair = [localize(events["dry"], events["released"], tr) for tr in (1, 2)]
                    focus = []
                    for r in pair:
                        # Trace both the earliest changed derivative call and
                        # the first changed returned local geometry, if any.
                        probes = []
                        for suffix in ("_path", "_coordinates"):
                            changed = [(r[phase+suffix], i, phase) for i, phase in enumerate(("hunt","second"))
                                       if phase+suffix in r]
                            step, _, call = min(changed) if changed else (0, 0, "")
                            r["probe"+suffix] = call
                            probes.append(step)
                        focus.append(probes)
            for tr, (path_step, geom_step) in enumerate(focus, start=1):
                summary = dict(drug=drug, seed=seed, trajectory=tr, **pair[tr-1])
                for s, suffix in ((path_step, "_path"), (geom_step,"_coordinates")):
                    if not s:
                        continue
                    evidence = [read_evals(raw/f"{drug}_{seed}_{arm}_focus.trace", tr, s, summary["probe"+suffix])
                                for arm in ("dry", "released")]
                    details = summarize_evals(*evidence)
                    if suffix == "_path":
                        summary.update(details)
                    else:
                        summary["first_changed_local_geometry_step"] = s
                        summary["first_changed_local_geometry_rms_A"] = details["local_end_same_atom_rms_A"]
                cases.append(summary)
            rows(out/"runs.csv", cells)
            rows(out/"divergences.csv", cases)
    meta = {"source": "ae19b9fc52e7501d356e06cdc3eed296021d9147", "observer_sha256": sha(binary),
            "receptor_sha256": {k: sha(v) for k, v in paths.items()},
            "ligand_sha256": {k: sha(work/"prepared"/(k+".pdbqt")) for k in ("at9283","danusertib")}}
    (out/"inputs.json").write_text(json.dumps(meta, indent=2)+"\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("work", "far_receptor", "observer", "output"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    run(*(p.resolve() for p in (a.work, a.far_receptor, a.observer, a.output)))
