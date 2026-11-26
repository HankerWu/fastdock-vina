"""Replay matched global searches for each measured kinase inhibitor."""

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess

STATES = {"CHEK1": ("2E9P", "2YWP"), "AURKA": ("2X6E", "2J4Z")}
SEEDS = (60913, 60943, 61211, 61231)
ARMS = {
    "native": [],
    "balanced": ["--balanced_torsion_search"],
    "prepared": ["--prepared_torsion_starts"],
}


def run(work, binary, arms, jobs):
    data = work / "prepared"
    rows = list(csv.DictReader((work / "panel.csv").open()))
    out = work / "pools"
    out.mkdir(exist_ok=True)

    def dock(case):
        row, state, seed, arm = case
        stem = f"{state}_{row['cid']}_{seed}_{arm}"
        dest = out / (stem + ".pdbqt")
        center = json.loads((data / (state + ".json")).read_text())["center"]
        box = []
        for ax, pos in zip("xyz", center):
            box.extend(["--center_" + ax, str(pos), "--size_" + ax, "26"])
        with (out / (stem + ".log")).open("w") as f:
            subprocess.run(
                [
                    binary,
                    "--receptor",
                    data / (state + ".pdbqt"),
                    "--ligand",
                    data / (row["cid"] + ".pdbqt"),
                    "--pose_selection",
                    "late",
                    "--cpu",
                    "1",
                    "--exhaustiveness",
                    "2",
                    "--max_evals",
                    "25000",
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--search_counters",
                    "--seed",
                    str(seed),
                    "--out",
                    dest,
                    *box,
                    *ARMS[arm],
                ],
                check=True,
                stdout=f,
                stderr=subprocess.STDOUT,
            )

    cases = [
        (r, s, k, a)
        for k in SEEDS
        for r in rows
        for s in STATES[r["target"]]
        for a in arms
    ]
    with ThreadPoolExecutor(jobs) as e:
        list(e.map(dock, cases))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--arms", nargs="+", choices=ARMS, default=list(ARMS))
    p.add_argument("--jobs", type=int, default=3)
    a = p.parse_args()
    run(a.work.resolve(), a.binary.resolve(), a.arms, a.jobs)
