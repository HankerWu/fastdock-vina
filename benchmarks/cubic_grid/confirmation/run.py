"""Independent CPU-cost calibration, reciprocal docking and official Vina comparison."""
import argparse
import csv
import itertools
import json
import os
from pathlib import Path
import statistics
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.spatial.transform import Rotation

SEEDS = (68041, 95087)
CALIBRATION_SEED = 90709
BASE_CAP = 24000
CORE_TARGETS = ("EPAS1", "FKBP1A", "RBP4", "FABP4", "DYRK1A")


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def experiment(job):
    fixtures, out, binary, query, receiver, arm, seed, cap = job
    key = f"{query['pdb']}-{receiver['pdb']}-{arm}-{seed}-{cap}"
    pose = out / (key + ".pdbqt")
    args = [str(binary), "--receptor", str(fixtures / receiver["pdb"] / "receptor.pdbqt"),
            "--ligand", str(fixtures / query["pdb"] / "ligand.pdbqt"),
            "--out", str(pose), "--cpu", "1", "--exhaustiveness", "3",
            "--max_evals", str(cap), "--num_modes", "9", "--energy_range", "100", "--seed", str(seed)]
    for axis, center in zip("xyz", receiver["center"]):
        args.extend(("--size_" + axis, "24", "--center_" + axis, str(center)))
    if arm != "official":
        args.extend(("--pose_selection", "late", "--refine_candidates", "15", "--search_counters"))
    if arm == "cubic":
        args.append("--cubic_search_grid")
    args.extend(("--spacing", "0.25" if arm == "dense" else "0.375"))
    t0, cpu0 = time.monotonic(), os.times()
    with (out / (key + ".log")).open("w") as log:
        subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, check=True)
    elapsed, cpu1 = time.monotonic()-t0, os.times()
    cpu = (cpu1.children_user-cpu0.children_user) + (cpu1.children_system-cpu0.children_system)
    return dict(target=query["target"], query=query["pdb"], receiver=receiver["pdb"], arm=arm,
                seed=seed, cap=cap, cpu_seconds=round(cpu, 3), wall_seconds=round(elapsed, 4),
                pool=pose.name)


def reference(query, receiver):
    # This path uses scipy's quaternion fit, independently of the primary SVD.
    keys = sorted(query["ca"].keys() & receiver["ca"].keys())
    keys = [k for k in keys if all(np.linalg.norm(np.asarray(e["ca"][k])-e["center"]) < 14
                                   for e in (query, receiver))]
    if len(keys) < 12:
        raise ValueError("Insufficient matched pocket C-alpha atoms")
    x, y = [np.asarray([e["ca"][key] for key in keys]) for e in (query, receiver)]
    rot, _ = Rotation.align_vectors(y-y.mean(0), x-x.mean(0))
    original = np.asarray(query["mapped"])
    return rot.apply((original-x.mean(0)).reshape(-1,3)).reshape(original.shape) + y.mean(0)


def measure(path, serial, expected):
    rmsd, coords = [], {}
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            coords[int(line[6:11])] = [float(line[start:start+8]) for start in (30,38,46)]
        elif line.startswith("ENDMDL"):
            observed = np.asarray([coords[i] for i in serial])
            distances = np.linalg.norm(observed[None, :, :] - expected, axis=2)
            rmsd.append(float(np.sqrt(np.mean(distances**2, axis=1)).min()))
            coords.clear()
    if not rmsd:
        raise ValueError(f"Empty pose pool {path}")
    return rmsd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=("calibrate", "confirm"))
    ap.add_argument("fixtures", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("vina", type=Path)
    ap.add_argument("--official", type=Path)
    ap.add_argument("--jobs", type=int, default=1)
    args = ap.parse_args()
    if args.stage == "confirm" and args.official is None:
        ap.error("confirm requires --official")
    out, fixtures = args.out.resolve(), args.fixtures.resolve()
    out.mkdir(parents=True, exist_ok=True)
    entries = [json.loads(p.read_text()) for p in sorted(fixtures.glob("*/preparation.json"))]
    entries = {e["pdb"]: e for e in entries if e["status"] == "ready"}
    jobs = []
    for target in CORE_TARGETS:
        pair = [e for e in entries.values() if e["target"] == target]
        if len(pair) != 2:
            raise ValueError(f"Expected two complete structures for {target}")
        if args.stage == "calibrate":
            directions, seeds = [tuple(pair)], [CALIBRATION_SEED]
            arms = [("native", BASE_CAP), ("cubic", BASE_CAP), ("native", 2*BASE_CAP)]
        else:
            directions, seeds = list(itertools.permutations(pair, 2)), SEEDS
            matched_cap = json.loads((out / "calibration.json").read_text())["native_matched_cap"]
            arms = [("native", BASE_CAP), ("cubic", BASE_CAP), ("dense", BASE_CAP),
                    ("native_matched", matched_cap), ("official", BASE_CAP)]
        for (query, receiver), seed, (arm, cap) in itertools.product(directions, seeds, arms):
            binary = args.official if arm == "official" else args.vina
            jobs.append((fixtures, out, binary.resolve(), query, receiver, arm, seed, cap))
    # Interleave arms, then reverse alternate case groups to balance warm cache effects.
    width = len(arms)
    jobs = [j for i in range(0,len(jobs),width)
            for j in (jobs[i:i+width] if i//width % 2 == 0 else reversed(jobs[i:i+width]))]
    with ProcessPoolExecutor(args.jobs) as executor:
        rows = list(executor.map(experiment, jobs))
    write_csv(out / (args.stage+"_timings.csv"), rows)
    if args.stage == "calibrate":
        totals = [sum(r["cpu_seconds"] for r in rows if r["arm"] == arm and r["cap"] == cap)
                  for arm,cap in arms]
        slope = (totals[2]-totals[0])/BASE_CAP
        if slope <= 0:
            raise ValueError("Calibration did not resolve the incremental native cost")
        cap = max(1000, 1000*round((BASE_CAP+(totals[1]-totals[0])/slope)/1000))
        (out/"calibration.json").write_text(json.dumps(dict(native_24k_cpu=totals[0],
            cubic_24k_cpu=totals[1], native_48k_cpu=totals[2], native_matched_cap=cap),indent=2)+"\n")
        return
    for row in rows:
        q, rec = entries[row["query"]], entries[row["receiver"]]
        values = measure(out/row["pool"], q["serial"], reference(q, rec))
        row.update(first=values[0], best5=min(values[:5]), best=min(values), modes=len(values))
    write_csv(out/"confirmation.csv",rows)
    for arm, _ in arms:
        group=[r for r in rows if r["arm"]==arm]
        print(arm,len(group),[sum(r[k]<=2 for r in group) for k in ("first","best5","best")],
              round(sum(r["cpu_seconds"] for r in group),2),round(statistics.median(r["first"] for r in group),3))


if __name__ == "__main__":
    main()
