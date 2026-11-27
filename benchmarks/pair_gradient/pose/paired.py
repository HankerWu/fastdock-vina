"""Paired one-thread redocking and cross-docking, with fixed structures and seeds."""

import argparse
import csv
import itertools
import json
import math
import random
import re
import subprocess
import time
from pathlib import Path

import gemmi
import numpy as np
from prepare import CASES, base
from scipy.spatial import cKDTree

ARMS = ("upstream", "parent", "linear", "hermite")
SEEDS = (75467, 98251)
CAPS = (3000, 12000)
RADII = {
    "C": 1.70,
    "A": 1.70,
    "N": 1.55,
    "NA": 1.55,
    "OA": 1.52,
    "O": 1.52,
    "S": 1.80,
    "SA": 1.80,
    "P": 1.80,
    "F": 1.47,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
}


def atoms(lines):
    data = [
        (int(s[6:11]), base.coords(s), s.split()[-1])
        for s in lines
        if s.startswith(("ATOM  ", "HETATM")) and s.split()[-1] not in {"H", "HD"}
    ]
    return {serial: (xyz, RADII[typ]) for serial, xyz, typ in data}


def renin_correspondence(work, meta):
    # Depositors used two incompatible author-numbering systems. Align observed
    # polymer sequences to one common deposited full-length mature REN sequence.
    ref = gemmi.read_structure(str(work / "raw/1RNE.pdb"))
    ref.setup_entities()
    seq = ref.entities[0].full_sequence
    for pdb, item in meta.items():
        if item["target"] != "REN":
            continue
        structure = gemmi.read_structure(str(work / "raw" / (pdb + ".pdb")))
        structure.setup_entities()
        poly = structure[0][item["chain"]].get_polymer()
        alignment = gemmi.align_sequence_to_polymer(
            seq, poly, gemmi.PolymerType.PeptideL
        )
        common = {}
        i = j = 0
        for count, op in re.findall(r"(\d+)([MID])", alignment.cigar_str()):
            for _ in range(int(count)):
                if op == "M" and seq[i] == poly[j].name:
                    ca = poly[j].find_atom("CA", "*")
                    if ca is not None:
                        common[f"{i + 1}:{poly[j].name}"] = list(ca.pos)
                i += op in "MI"
                j += op in "MD"
        if i != len(seq) or j != len(poly) or alignment.calculate_identity() < 99:
            raise ValueError("Unreliable REN sequence correspondence for " + pdb)
        item["ca"] = common


def inspect(directory, query, receiver, row):
    rot, trans, ca_rmsd, ca_count = base.alignment(query, receiver)
    refs = np.array(query["mapped"]) @ rot + trans
    rec = atoms(
        (directory.parent / receiver["pdb"] / "receptor.pdbqt").read_text().splitlines()
    )
    rxyz = np.array([v[0] for v in rec.values()])
    rr = np.array([v[1] for v in rec.values()])
    tree = cKDTree(rxyz)
    poses = []
    for idx, block in enumerate(
        (directory / (row["key"] + ".pdbqt")).read_text().split("ENDMDL")[:-1], 1
    ):
        lines = block.splitlines()
        heavy = atoms(lines)
        xyz = np.array([heavy[i][0] for i in query["serial"]])
        rad = np.array([heavy[i][1] for i in query["serial"]])
        rms = np.sqrt(np.square(refs - xyz).sum(axis=2).mean(axis=1)).min()
        clashes = 0
        for i, neighbors in enumerate(tree.query_ball_point(xyz, 3.0)):
            d = np.linalg.norm(xyz[i] - rxyz[neighbors], axis=1)
            clashes += np.count_nonzero(d < 0.75 * (rr[neighbors] + rad[i]))
        poses.append(
            {
                "key": row["key"],
                "mode": idx,
                "rmsd_A": rms,
                "clash_pairs": clashes,
                "minimum_heavy_distance_A": float(tree.query(xyz)[0].min()),
                "score": float(
                    next(
                        l.split()[3]
                        for l in lines
                        if l.startswith("REMARK VINA RESULT:")
                    )
                ),
            }
        )
    if not poses:
        raise ValueError(row["key"])
    row.update(
        ca_rmsd_A=ca_rmsd,
        ca_count=ca_count,
        modes=len(poses),
        top1_A=poses[0]["rmsd_A"],
        best5_A=min(p["rmsd_A"] for p in poses[:5]),
        oracle_A=min(p["rmsd_A"] for p in poses),
        top1_clash_pairs=poses[0]["clash_pairs"],
    )
    return poses


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def counts(log):
    stats = {}
    for line in log.splitlines():
        if "Search trajectory:" in line or "Search counters:" in line:
            stats.update((k, int(v)) for k, v in re.findall(r"(\w+)=(\d+)", line))
    keys = (
        "hunt_calls",
        "hunt_evals",
        "second_calls",
        "second_evals",
        "explicit_inputs",
        "explicit_calls",
        "explicit_evals",
    )
    return {k: stats.get(k, "") for k in keys}


def execute(work, meta, s, r, seed, cap, arm, binary, resume, replay=None):
    directory = work / "runs"
    key = f"{s}_{r}_{seed}_{cap}_{arm}"
    out = directory / (key + ".pdbqt")
    log = directory / (key + ".log")
    wall = directory / (key + ".time")
    common = [
        "--receptor",
        str(work / r / "receptor.pdbqt"),
        "--ligand",
        str(work / s / "ligand.pdbqt"),
        "--cpu",
        "1",
        "--exhaustiveness",
        "1",
        "--max_evals",
        str(cap),
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--seed",
        str(seed),
    ]
    for axis, value in zip("xyz", meta[r]["center"]):
        common += ["--center_" + axis, str(value), "--size_" + axis, "24"]
    if not resume or not wall.exists():
        cmd = [str(binary.resolve())] + common + ["--out", str(out)]
        if arm not in ("upstream", "release"):
            cmd += ["--search_counters"]
        start = time.perf_counter()
        run = subprocess.run(cmd, capture_output=True, text=True, check=True)
        elapsed = time.perf_counter() - start
        log.write_text(run.stdout + run.stderr)
        wall.write_text(str(elapsed) + "\n")
    parity = ""
    count_log = log.read_text()
    if arm == "upstream" and replay:
        counts_path = directory / (key + ".counts")
        if not resume or not counts_path.exists():
            comparison = directory / (key + ".replay")
            run = subprocess.run(
                [str(replay.resolve())] + common + ["--out", str(comparison)],
                capture_output=True,
                text=True,
                check=True,
            )
            if out.read_bytes() != comparison.read_bytes():
                raise RuntimeError("upstream instrumented replay differs for " + key)
            counts_path.write_text(run.stdout + run.stderr)
            comparison.unlink()
        count_log = counts_path.read_text()
        parity = "identical"
    row = dict(
        key=key,
        target=meta[s]["target"],
        query=s,
        receiver=r,
        seed=seed,
        cap=cap,
        arm=arm,
        kind="cognate" if s == r else "cross",
        seconds=float(wall.read_text()),
        upstream_count_replay=parity,
        **counts(count_log),
    )
    poses = inspect(directory, meta[s], meta[r], row)
    print(row["key"], round(row["seconds"], 3), flush=True)
    return row, poses


def calibrate(cells, query, receiver, arm):
    timing = {
        (r["arm"], int(r["cap"])): float(r["seconds"])
        for r in cells
        if r["query"] == query
        and r["receiver"] == receiver
        and int(r["seed"]) == SEEDS[0]
    }
    low, high = CAPS
    target = timing["parent", high]
    slow, fast = timing[arm, high], timing[arm, low]
    if slow <= target:
        return high
    if fast >= target or fast >= slow:
        return low
    return math.floor(low + (high - low) * (target - fast) / (slow - fast))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    for arm in (*ARMS, "shift", "release"):
        p.add_argument("--" + arm, type=Path)
    p.add_argument(
        "--counts-replay",
        type=Path,
        help="same upstream source with passive eval counters",
    )
    p.add_argument("--only-target", choices=("ADORA2A", "REN"))
    p.add_argument(
        "--resume",
        action="store_true",
        help="read already completed run logs, poses and timings",
    )
    p.add_argument(
        "--phase", choices=("primary", "shift", "matched", "release"), default="primary"
    )
    args = p.parse_args()
    work = args.work.resolve()
    (work / "runs").mkdir(exist_ok=True)
    meta = {
        pdb: json.loads((work / pdb / "preparation.json").read_text()) for pdb in CASES
    }
    renin_correspondence(work, meta)
    pairs = [
        (s, r)
        for s, r in itertools.product(CASES, repeat=2)
        if meta[s]["target"] == meta[r]["target"]
        and meta[r]["receptor_status"] == "ready"
        and args.only_target in (None, meta[s]["target"])
    ]
    tasks = []
    if args.phase == "matched":
        cells = []
        for filename in ("cells.csv", "shift_cells.csv"):
            with (work / filename).open() as f:
                cells.extend(csv.DictReader(f))
        for s, r in pairs:
            for arm in ("linear", "hermite", "shift"):
                tasks.append((s, r, SEEDS[1], calibrate(cells, s, r, arm), arm))
        random.Random(20261127).shuffle(tasks)
    else:
        arms = ARMS if args.phase == "primary" else (args.phase,)
        for target in dict.fromkeys(meta[s]["target"] for s, _ in pairs):
            cells = [
                (s, r, seed, cap)
                for s, r in pairs
                if meta[s]["target"] == target
                for seed, cap in itertools.product(SEEDS, CAPS)
            ]
            random.Random(20261127).shuffle(cells)
            for i, cell in enumerate(cells):
                k = i % len(arms)
                tasks.extend((*cell, arm) for arm in arms[k:] + arms[:k])
    cell_rows, modes = [], []
    for s, r, seed, cap, arm in tasks:
        binary = getattr(args, arm)
        if binary is None:
            p.error("missing binary path for " + arm)
        reuse = args.resume or (args.phase == "matched" and cap in CAPS)
        row, poses = execute(
            work, meta, s, r, seed, cap, arm, binary, reuse, args.counts_replay
        )
        cell_rows.append(row)
        modes.extend(poses)
    prefix = "" if args.phase == "primary" else args.phase + "_"
    write(work / (prefix + "cells.csv"), cell_rows)
    write(work / (prefix + "modes.csv"), modes)


if __name__ == "__main__":
    main()
