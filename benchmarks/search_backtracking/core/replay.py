"""Run one hunt and save from each frozen start on each arm, with stage coordinates."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "geometry_comparison"
ARMS = ("native", "halving", "quadratic")


def pose(path, serial):
    xyz = {int(line[6:11]): [float(line[start:start + 8]) for start in (30, 38, 46)]
           for line in path.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))}
    return np.array([xyz[i] for i in serial])


def write_rows(path, rows):
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def rejected_profiles(log, row):
    matches = re.findall(r"Exhausted line (\d+) f0=(\S+) pg=(\S+) delta=(\S+) last_slope=(\S+) probes=(\S+)", log)
    if not matches:
        return []
    stages = ["hunt"] * row["hunt_exhausted"] + ["save"] * row["save_exhausted"]
    if len(stages) != len(matches):
        raise ValueError("incomplete exhaustion profile")
    profiles = []
    for stage, (index, f0, pg, delta, slope, extra) in zip(stages, matches):
        profile = {k: row[k] for k in ("pdb", "target", "start")}
        profile.update(stage=stage, line=int(index), f0=f0, initial_slope=pg, last_slope=slope)
        # delta_9 spans the last native trial at alpha = 2^-9. The remaining
        # ten values come from diagnostic evaluations at 2^-10 ... 2^-19.
        for exponent, value in enumerate([delta, *extra.split(",")], start=9):
            profile[f"delta_{exponent}"] = value
        profiles.append(profile)
    return profiles


def run(binary, work, arms):
    work.mkdir(parents=True)
    for archive, dest in [("pose/inputs.tar.xz", "fixtures"), ("basins/starts_and_relaxations.tar.xz", "starts")]:
        with tarfile.open(SOURCE / archive) as stream:
            stream.extractall(work / dest, filter="data")
    raw = work / "raw"
    raw.mkdir()
    rows, profiles = [], []
    for directory in sorted((work / "starts").iterdir()):
        fixture = work / "fixtures" / directory.name
        info = json.loads((fixture / "metadata.json").read_text())
        refs = np.asarray(info["mapped"])
        def rms(coords):
            return np.sqrt(np.sum((coords - refs) ** 2, axis=2).mean(axis=1)).min()
        for label in ("projected", "63419", "63421"):
            start = directory / (label + ".pdbqt")
            input_coords = pose(start, info["serial"])
            for arm in arms:
                prefix = raw / f"{directory.name}-{label}-{arm}"
                args = [binary, "--receptor", fixture / "receptor.pdbqt", "--ligand", start,
                        "--cpu", "1", "--seed", "63419", "--exhaustiveness", "1", "--max_evals", "1",
                        "--num_modes", "1", "--search_counters", "--out", str(prefix) + ".final.pdbqt"]
                for axis, center in zip("xyz", info["center"]):
                    args += ["--center_" + axis, str(center), "--size_" + axis, "24"]
                log = subprocess.run(args, check=True, capture_output=True, text=True,
                                     env={**os.environ, "FASTDOCK_RESEARCH_BACKTRACK": arm,
                                          "FASTDOCK_FIRST_HUNT_OUT": str(prefix)}).stdout
                (raw / (prefix.name + ".log")).write_text(log)
                counters = re.search(r"hunt_calls=1 hunt_evals=(\d+) second_calls=1 second_evals=(\d+)", log)
                if not counters:
                    raise ValueError(f"expected exactly one hunt and save in {prefix}")
                row = dict(pdb=directory.name, target=info["target"], start=label, arm=arm)
                for stage in ("start", "hunt", "save"):
                    coords = pose(raw / (prefix.name + "." + stage + ".pdbqt"), info["serial"])
                    if stage == "start" and not np.array_equal(coords, input_coords):
                        raise ValueError("Replayed coordinates differ from input: " + str(prefix))
                    row[stage + "_rmsd_A"] = round(rms(coords), 6)
                for stage, evaluations in zip(("hunt", "save"), counters.groups()):
                    match = re.search(r"Backtracking: stage=" + stage + r" searches=(\d+) rejected=(\d+) exhausted=(\d+)", log)
                    searches, rejected, exhausted = map(int, match.groups())
                    if int(evaluations) != 1 + searches + rejected - exhausted:
                        raise ValueError("trial count mismatch: " + str(prefix))
                    for key, value in zip(("evals", "searches", "rejected", "exhausted"),
                                          (int(evaluations), searches, rejected, exhausted)):
                        row[stage + "_" + key] = value
                rows.append(row)
                profiles.extend(rejected_profiles(log, row))
    if len(rows) != 24 * len(arms):
        raise ValueError("Incomplete eight-target fixed-start panel")
    write_rows(work / "rows.csv", rows)
    if profiles:
        write_rows(work / "exhausted_lines.csv", profiles)
    with tarfile.open(work / "raw.tar.xz", "w:xz") as bundle:
        bundle.add(raw, arcname="raw")
    print("replay_binary_sha256", hashlib.sha256(binary.read_bytes()).hexdigest())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--only", choices=("native", "quarter"), help="run a single diagnostic control")
    a = p.parse_args()
    run(a.binary.resolve(), a.work.resolve(), (a.only,) if a.only else ARMS)
