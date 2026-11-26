"""Run the six matched alkylbenzene cross-docking queries against each receptor state."""

import argparse
import itertools
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ARMS = ("native", "single", "multi")
SEEDS = (161803, 141421)
STATES = ("A", "B", "C")


def run(binary, work, cell):
    state, code, seed, arm = cell
    prep = work / "prepared"
    stem = work / "runs" / f"{state}_{code}_{seed}_{arm}"
    output = stem.with_suffix(".pdbqt")
    command = [
        binary,
        "--ligand",
        prep / f"{code}_ligand.pdbqt",
        "--receptor",
        prep / f"{state}_receptor.pdbqt",
        "--config",
        prep / "box.txt",
        "--cpu",
        "1",
        "--seed",
        str(seed),
        "--exhaustiveness",
        "4",
        "--max_evals",
        "50000",
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--pose_selection",
        "late",
        "--search_counters",
        "--out",
        output,
    ]
    if arm != "native":
        command += ["--search_anchor", prep / f"{code}.{arm}"]
    with stem.with_suffix(".log").open("w") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--jobs", type=int, default=3)
    a = p.parse_args()
    (a.work / "runs").mkdir(exist_ok=True)
    queries = sorted(json.loads((a.work / "prepared/manifest.json").read_text()))
    with ThreadPoolExecutor(a.jobs) as pool:
        list(
            pool.map(
                lambda c: run(a.binary.resolve(), a.work, c),
                itertools.product(STATES, queries, SEEDS, ARMS),
            )
        )
