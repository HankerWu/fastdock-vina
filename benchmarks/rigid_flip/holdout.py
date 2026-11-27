"""Prepare the preregistered human structures and replay all four blind docking arms."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random
import re
import subprocess
import sys

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH / "search_backtracking/pose"))
import backtracking_panel as panel

paired = panel.paired
paired.RADII["Mg"] = 1.73  # Same validation radius as the existing bond_distance/pose panel.
CASES = {
    "1T48": ("PTPN1", "BB3", "A"),
    "1T49": ("PTPN1", "892", "A"),
    "3ERT": ("ESR1", "OHT", "A"),
    "1SJ0": ("ESR1", "E4D", "A"),
    "4UND": ("PARP1", "2YQ", "A"),
}
SEEDS = (202711, 202727)
CAPS = (3000, 12000)
ARMS = ("upstream", "parent", "always", "adaptive")


def args_for(work, meta, q, r, seed, cap, out):
    args = ["--receptor", str(work / r / "receptor.pdbqt"), "--ligand", str(work / q / "ligand.pdbqt"),
            "--cpu", "1", "--exhaustiveness", "1", "--seed", str(seed), "--max_evals", str(cap),
            "--num_modes", "9", "--energy_range", "100", "--out", str(out)]
    for k, c in zip("xyz", meta[r]["center"]):
        args.extend(["--center_" + k, str(c), "--size_" + k, "24"])
    return args


def run(work, binary, replay, meta, q, r, seed, cap, arm):
    prefix = work / "runs" / f"{q}_{r}_{seed}_{cap}_{arm}"
    out = prefix.with_suffix(".pdbqt")
    log = prefix.with_suffix(".log")
    stats = prefix.with_suffix(".json")
    if stats.exists():
        return
    cli = args_for(work, meta, q, r, seed, cap, out)
    cost = prefix.with_suffix(".cost")
    subprocess.run([sys.executable, "-S", str(Path(__file__).with_name("process_cost.py")), str(log),
                    str(cost), str(binary.resolve()), *cli, *([] if arm == "upstream" else ["--search_counters"])],
                   check=True)
    observed = log.read_text()
    row = dict(key=prefix.name, target=CASES[q][0], query=q, receiver=r, seed=seed, cap=cap, arm=arm,
               kind="cognate" if q == r else "cross", **json.loads(cost.read_text()), sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
               upstream_counter_parity="")
    if arm == "upstream":
        check = prefix.with_suffix(".replay")
        result = subprocess.run([str(replay.resolve()), *args_for(work, meta, q, r, seed, cap, check)],
                                capture_output=True, text=True, check=True)
        if check.read_bytes() != out.read_bytes():
            raise ValueError(f"Upstream counter replay changed poses for {prefix.name}")
        check.unlink()
        observed = result.stdout + result.stderr
        prefix.with_suffix(".counts").write_text(observed)
        row["upstream_counter_parity"] = "identical"
    row.update(paired.counts(observed))
    flip = re.search(r"\brigid_reorientations=(\d+)", observed)
    row["rigid_reorientations"] = int(flip[1]) if flip else (0 if arm == "parent" else "")
    if any(v == "" for v in paired.counts(observed).values()):
        raise ValueError(f"Missing actual eval count for {prefix.name}")
    stats.write_text(json.dumps(row, indent=2) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=["prepare", "dock", "analyze"])
    p.add_argument("work", type=Path)
    p.add_argument("--ids", nargs="+", choices=CASES)
    for name in (*ARMS, "counts"):
        p.add_argument("--" + name, type=Path)
    p.add_argument("--study", choices=("primary", "reduced", "newseeds"), default="primary")
    a = p.parse_args()
    work = a.work.resolve()
    if a.action == "prepare":
        panel.CASES = CASES
        for code in a.ids or CASES:
            panel.prepare_input(work, code)
        return
    meta = {s: json.loads((work / s / "preparation.json").read_text()) for s in CASES}
    ready = {s for s in CASES if meta[s]["receptor_status"] == "ready"}
    seeds = (91373, 91411) if a.study == "newseeds" else SEEDS
    caps = {"primary": CAPS, "reduced": (2550, 10200), "newseeds": (12000,)}[a.study]
    cells = []
    excluded = []
    for q, r in itertools.product(sorted(CASES), sorted(ready)):
        if CASES[q][0] != CASES[r][0]:
            continue
        try:
            paired.base.alignment(meta[q], meta[r])
        except ValueError as err:
            excluded.append(dict(q=q, r=r, reason=str(err)))
            continue
        cells.extend((q, r, s, c) for s, c in itertools.product(seeds, caps))
    (work / "alignment_exclusions.json").write_text(json.dumps(excluded, indent=2) + "\n")
    (work / "runs").mkdir(exist_ok=True)
    if a.action == "dock":
        order = random.Random(66759)
        order.shuffle(cells)
        for cell in cells:
            arms = [name for name in ARMS if getattr(a, name)]
            order.shuffle(arms)
            for arm in arms:
                run(work, getattr(a, arm), a.counts, meta, *cell, arm)
        return
    rows, poses = [], []
    for path in sorted((work / "runs").glob("*.json")):
        row = json.loads(path.read_text())
        for k in ("hunt_evals", "second_evals", "explicit_evals"):
            if not isinstance(row[k], int):
                raise ValueError((path, k))
        models = paired.inspect(work / "runs", meta[row["query"]], meta[row["receiver"]], row)
        if not models:
            raise ValueError(path)
        rows.append(row)
        poses.extend(models)
    paired.write(work / "cells.csv", rows)
    paired.write(work / "modes.csv", poses)


if __name__ == "__main__":
    main()
