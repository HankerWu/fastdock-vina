"""Interleaved, paired docking and upstream counter replay on one pinned core."""

import argparse
import hashlib
import itertools
import json
import os
import random
import re
import shutil
from pathlib import Path

import numpy as np
from backtracking_panel import ARMS, CAPS, CASES, SEEDS, paired, write_sources
from scipy.spatial.distance import cdist


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geometry(work, meta):
    rows = []
    for q, r in itertools.product(CASES, repeat=2):
        a, b = meta[q], meta[r]
        if a["target"] != b["target"] or b["receptor_status"] != "ready":
            continue
        rot, shift, rms, count = paired.base.alignment(a, b)
        receptor = paired.atoms((work / r / "receptor.pdbqt").read_text().splitlines())
        ligand = paired.atoms((work / q / "ligand.pdbqt").read_text().splitlines())
        d = cdist(
            np.array(a["mapped"][0]) @ rot + shift, [v[0] for v in receptor.values()]
        )
        rad = np.array([ligand[i][1] for i in a["serial"]])[:, None] + [
            v[1] for v in receptor.values()
        ]
        rows.append(
            {
                "query": q,
                "receiver": r,
                "pocket_ca": count,
                "alignment_rmsd_A": rms,
                "reference_clashes": int((d < 0.75 * rad).sum()),
                "reference_min_distance_A": d.min(),
            }
        )
    paired.write(work / "reference_geometry.csv", rows)
    return [(row["query"], row["receiver"]) for row in rows]


def run(args):
    work = args.work.resolve()
    os.sched_setaffinity(0, {args.cpu})
    meta = {
        code: json.loads((work / code / "preparation.json").read_text())
        for code in CASES
    }
    write_sources(work)
    pairs = geometry(work, meta)
    cells = [(q, r, s, c) for q, r in pairs for s, c in itertools.product(SEEDS, CAPS)]
    random.Random(SEEDS[0]).shuffle(cells)
    (work / "runs").mkdir(exist_ok=True)
    bins = {arm: getattr(args, arm).resolve() for arm in ARMS}
    rows, modes, verify = [], [], []
    for i, (q, r, s, c) in enumerate(cells):
        start = i % len(ARMS)
        for arm in ARMS[start:] + ARMS[:start]:
            row, poses = paired.execute(
                work, meta, q, r, s, c, arm, bins[arm], args.resume, args.counts
            )
            log = (work / "runs" / (row["key"] + ".log")).read_text()
            for stage in ("hunt", "save"):
                found = re.search(
                    r"Backtracking: stage="
                    + stage
                    + r" searches=(\d+) rejected=(\d+) exhausted=(\d+)",
                    log,
                )
                for j, metric in enumerate(("searches", "rejected", "exhausted"), 1):
                    row[stage + "_" + metric] = int(found[j]) if found else ""
                if arm != "upstream" and not found:
                    raise ValueError("Missing stage counters: " + row["key"])
            for k in ("hunt_evals", "second_evals", "explicit_evals"):
                if row[k] == "":
                    raise ValueError("Missing work counter: " + row["key"])
            rows.append(row)
            modes.extend(poses)
    for target in sorted({a["target"] for a in meta.values()}):
        q = next(q for q, r in pairs if q == r and meta[q]["target"] == target)
        for arm in ("parent", "quadratic"):
            s, c = SEEDS[0], CAPS[1]
            row = next(
                v
                for v in rows
                if (v["query"], v["receiver"], v["seed"], v["cap"], v["arm"])
                == (q, q, s, c, arm)
            )
            dest = work / "replay"
            dest.mkdir(exist_ok=True)
            (dest / "runs").mkdir(exist_ok=True)
            if not (dest / q).exists():
                shutil.copytree(work / q, dest / q)
            repeat, _ = paired.execute(dest, meta, q, q, s, c, arm, bins[arm], False)
            same = digest(dest / "runs" / (row["key"] + ".pdbqt")) == digest(
                work / "runs" / (row["key"] + ".pdbqt")
            )
            if not same:
                raise ValueError("Output replay differed: " + row["key"])
            verify.append(
                {"key": row["key"], "identical": same, "seconds": repeat["seconds"]}
            )
    paired.write(work / "cells.csv", rows)
    paired.write(work / "modes.csv", modes)
    paired.write(work / "replays.csv", verify)
    info = {k: digest(v) for k, v in {**bins, "upstream_counters": args.counts}.items()}
    if args.shared_core:
        info["shared_core_binary"] = digest(args.shared_core)
    (work / "binaries.json").write_text(json.dumps(info, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    for name in (*ARMS, "counts"):
        p.add_argument("--" + name, required=True, type=Path)
    p.add_argument(
        "--shared-core",
        type=Path,
        help="compiled binary called by the three mode wrappers",
    )
    p.add_argument("--cpu", required=True, type=int)
    p.add_argument("--resume", action="store_true")
    run(p.parse_args())
