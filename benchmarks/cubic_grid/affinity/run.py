"""Paired equal-evaluation searches with identical native refinement and scoring."""

import argparse
import csv
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from time import monotonic

SEEDS = (76123, 76231)
STATES = {"SYK": "6ZCY", "FAK": "6YT6"}


def box_args(prepared, state):
    out = []
    for a, v in zip(
        "xyz", json.loads((prepared / f"{state}.json").read_text())["center"]
    ):
        out.extend(["--center_" + a, str(v), "--size_" + a, "26"])
    return out


def run(work, binary, jobs, arms, targets):
    prepared = work / "prepared"
    pools = work / "pools"
    pools.mkdir(exist_ok=True)
    panel = list(csv.DictReader((work / "panel.csv").open()))

    def dock(job):
        r, seed, arm = job
        state = STATES[r["target"]]
        stem = f"{state}_{r['cid']}_{seed}_{arm}"
        dest = pools / (stem + ".pdbqt")
        log = pools / (stem + ".log")
        if dest.exists() and dest.read_text().endswith("ENDMDL\n") and log.is_file():
            return
        tmp = dest.with_suffix(".working.pdbqt")
        cmd = [
            binary,
            "--receptor",
            prepared / (state + ".pdbqt"),
            "--ligand",
            prepared / (r["cid"] + ".pdbqt"),
            "--pose_selection",
            "late",
            "--cpu",
            "1",
            "--seed",
            str(seed),
            "--exhaustiveness",
            "2",
            "--max_evals",
            "12000",
            "--num_modes",
            "9",
            "--energy_range",
            "100",
            "--search_counters",
            "--out",
            tmp,
            *box_args(prepared, state),
        ]
        if arm == "cubic":
            cmd.append("--cubic_search_grid")
        t = monotonic()
        proc = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
        )
        log.write_text(proc.stdout + f"\nwall_s {monotonic() - t:.6f}\n")
        proc.check_returncode()
        assert tmp.read_text().endswith("ENDMDL\n")
        tmp.replace(dest)

    tasks = [
        (r, k, a) for k in SEEDS for r in panel if r["target"] in targets for a in arms
    ]
    with ThreadPoolExecutor(jobs) as ex:
        list(ex.map(dock, tasks))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument(
        "--arms", nargs="+", choices=("native", "cubic"), default=("native", "cubic")
    )
    p.add_argument("--targets", nargs="+", choices=STATES, default=STATES)
    a = p.parse_args()
    run(a.work.resolve(), a.binary.resolve(), a.jobs, a.arms, a.targets)
