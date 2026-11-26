"""Run equal-budget paired, full Late searches and retain every output and counter."""
import argparse
import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

STATES = {"BACE1": ("4DJW", "4GID"), "MAPK14": ("3GCP", "1KV2"),
          "ALK": ("7BTT",), "RET": ("7DUA",), "TNK2": ("7KP6",)}
SEEDS = (89513, 89527)
ARMS = ("armijo", "geometry")


def run(work, binary, jobs, cap, seeds):
    with (work / "panel.csv").open() as stream:
        panel = list(csv.DictReader(stream))
    out = work / "pools"
    out.mkdir(exist_ok=True)
    prepared = work / "prepared"

    def pair(job):
        state, row, seed = job
        arms = ARMS if seed == seeds[0] else ARMS[::-1]
        for arm in arms:
            name = f"{state}_{row['compound']}_{seed}_{arm}"
            path, log = out / (name + ".pdbqt"), out / (name + ".log")
            tmp = path.with_suffix(".tmp.pdbqt")
            command = [binary, "--receptor", prepared / (state + ".pdbqt"),
                       "--ligand", prepared / (row["compound"] + ".pdbqt"),
                       "--config", prepared / (state + ".box"), "--pose_selection", "late",
                       "--refinement_optimizer", arm, "--seed", str(seed), "--cpu", "1",
                       "--exhaustiveness", "2", "--max_evals", str(cap), "--num_modes", "9",
                       "--energy_range", "100", "--search_counters", "--out", tmp]
            with log.open("w") as stream:
                subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)
            if not tmp.read_bytes().endswith(b"ENDMDL\n"):
                raise ValueError(name)
            tmp.replace(path)

    tasks = [(s, r, k) for k in seeds for r in panel for s in STATES[r["target"]]]
    with ThreadPoolExecutor(jobs) as executor:
        list(executor.map(pair, tasks))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--cap", type=int, default=12000)
    parser.add_argument("--seeds", type=int, nargs="+", default=SEEDS)
    args = parser.parse_args()
    run(args.work.resolve(), args.binary.resolve(), args.jobs, args.cap, args.seeds)
