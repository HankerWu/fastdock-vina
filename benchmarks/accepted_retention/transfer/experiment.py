"""One-trajectory parent/upstream/one-slot docking on a frozen structural panel."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
import re
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

from inputs import CAPS, CASES, SEEDS, base, reference

ARMS = ("upstream", "parent", "gated", "unconditional")
KEYS = (
    "skipped",
    "checks",
    "records",
    "eligible",
    "replaced",
    "slots",
    "slot",
    "pending_energy",
    "previous_energy",
)


def inputs(work):
    meta = {p: json.loads((work / p / "preparation.json").read_text()) for p in CASES}
    geometries = []
    for q, r in itertools.product(CASES, repeat=2):
        a, b = meta[q], meta[r]
        if a["target"] != b["target"]:
            continue
        rot, t, rms, count = base.alignment(a, b)
        ra = reference.atoms((work / r / "receptor.pdbqt").read_text().splitlines())
        qa = reference.atoms((work / q / "ligand.pdbqt").read_text().splitlines())
        d = cdist(
            np.array(a["mapped"][0]) @ rot + t, np.array([v[0] for v in ra.values()])
        )
        rad = np.array([qa[i][1] for i in a["serial"]])[:, None] + np.array(
            [v[1] for v in ra.values()]
        )
        geometries.append(
            dict(
                query=q,
                receiver=r,
                aligned_ca=count,
                alignment_rmsd_A=rms,
                reference_clashes=int(np.sum(d < 0.75 * rad)),
                reference_min_distance_A=d.min(),
            )
        )
    reference.write(work / "reference_geometry.csv", geometries)
    return meta, [(r["query"], r["receiver"]) for r in geometries]


def run(work, meta, q, r, seed, cap, arm, bins, counts, reuse):
    row, poses = reference.execute(
        work, meta, q, r, seed, cap, arm, bins[arm], reuse, counts
    )
    vals = {}
    if arm in ("gated", "unconditional"):
        line = next(
            l
            for l in (work / "runs" / (row["key"] + ".log")).read_text().splitlines()
            if "Deferred hunt:" in l
        )
        vals = dict(re.findall(r"(\w+)=([-+.\w]+)", line))
        if set(KEYS) - set(vals):
            raise ValueError(vals)
    row.update({k: vals.get(k, "") for k in KEYS})
    return row, poses


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    for arm in (*ARMS, "counts", "release"):
        p.add_argument("--" + arm, type=Path)
    p.add_argument("--cpu", type=int, default=8)
    p.add_argument("--phase", choices=("primary", "release", "time"), default="primary")
    p.add_argument("--resume", action="store_true")
    a = p.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    work = a.work.resolve()
    (work / "runs").mkdir(exist_ok=True)
    meta, pairs = inputs(work)
    bins = {k: getattr(a, k) for k in (*ARMS, "release")}
    cells = [
        (*sr, seed, cap) for sr in pairs for seed, cap in itertools.product(SEEDS, CAPS)
    ]
    random.Random(214087).shuffle(cells)
    rows = []
    modes = []
    if a.phase in ("primary", "release"):
        arms = ARMS if a.phase == "primary" else ("release",)
        for i, (q, r, seed, cap) in enumerate(cells):
            j = i % len(arms)
            for arm in arms[j:] + arms[:j]:
                row, poses = run(
                    work, meta, q, r, seed, cap, arm, bins, a.counts, a.resume
                )
                if arm == "release":
                    peer = (
                        work
                        / "runs"
                        / (row["key"].replace("_release", "_upstream") + ".pdbqt")
                    )
                    row["identical_upstream"] = (
                        peer.read_bytes()
                        == (work / "runs" / (row["key"] + ".pdbqt")).read_bytes()
                    )
                rows.append(row)
                modes.extend(poses)
    else:
        with (work / "cells.csv").open() as f:
            raw = list(csv.DictReader(f))
        indexed = {
            (r["query"], r["receiver"], int(r["seed"]), int(r["cap"]), r["arm"]): r
            for r in raw
        }

        def total(arm, cost):
            return sum(sum(float(r[k]) for k in cost) for r in raw if r["arm"] == arm)

        for arm in ("gated", "unconditional"):
            excess = max(
                total(arm, cols) / total("parent", cols)
                for cols in [
                    ["seconds"],
                    ["hunt_evals", "second_evals", "explicit_evals"],
                ]
            )
            if excess <= 1.05:
                continue
            for q, r in pairs:

                def at(cap, method):
                    return float(indexed[q, r, SEEDS[0], cap, method]["seconds"])

                lo, hi = CAPS
                t0, t1, target = at(lo, arm), at(hi, arm), at(hi, "parent")
                cap = (
                    hi
                    if t1 <= target
                    else lo
                    if t0 >= target or t1 <= t0
                    else int(lo + (hi - lo) * (target - t0) / (t1 - t0))
                )
                row, poses = run(
                    work,
                    meta,
                    q,
                    r,
                    SEEDS[1],
                    cap,
                    arm,
                    bins,
                    None,
                    cap in CAPS or a.resume,
                )
                row.update(
                    parent_cap=hi,
                    calibration_low_s=t0,
                    calibration_high_s=t1,
                    calibration_parent_s=target,
                )
                rows.append(row)
                modes.extend(poses)
    prefix = "" if a.phase == "primary" else a.phase + "_"
    if rows:
        reference.write(work / (prefix + "cells.csv"), rows)
        reference.write(work / (prefix + "modes.csv"), modes)
    info = {
        arm: hashlib.sha256(path.read_bytes()).hexdigest()
        for arm, path in {**bins, "counts": a.counts}.items()
        if path
    }
    (work / (prefix + "binaries.json")).write_text(json.dumps(info, indent=2) + "\n")


if __name__ == "__main__":
    main()
