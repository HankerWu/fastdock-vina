"""Replay the frozen AD4 structural panel; all RMSDs use receptor coordinates."""

import argparse
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import random
import re
import subprocess
import time

import numpy as np
from prepare import ENTRIES, preparation


def write(path, rows):
    with path.open("w", newline="") as stream:
        out = csv.DictWriter(stream, rows[0], lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def run(args, job):
    query, receiver, seed, cap, arm = job
    binary = args.official if arm == "vina" else getattr(args, arm)
    prefix = args.out / "runs" / f"{query}-{receiver}-{seed}-{cap}-{arm}"
    cmd = [
        str(binary),
        "--ligand",
        str(args.work / query / "ligand.pdbqt"),
        "--scoring",
        "ad4",
        "--maps",
        str(args.work / "maps" / receiver / "receptor"),
        "--cpu",
        "1",
        "--exhaustiveness",
        "1",
        "--num_modes",
        "5",
        "--energy_range",
        "100",
        "--max_evals",
        str(cap),
        "--seed",
        str(seed),
        "--out",
        str(prefix.with_suffix(".pdbqt")),
    ]
    if arm == "vina":
        idx = cmd.index("--scoring")
        rec = json.loads((args.work / receiver / "preparation.json").read_text())
        cmd[idx : idx + 4] = [
            "--scoring",
            "vina",
            "--receptor",
            str(args.work / receiver / "receptor.pdbqt"),
        ]
        for axis, coord in zip("xyz", rec["center"]):
            cmd.extend(["--center_" + axis, f"{coord:.3f}", "--size_" + axis, "24"])
    if arm not in ("official", "vina"):
        cmd.append("--search_counters")
    with prefix.with_suffix(".log").open("w") as log:
        t = time.perf_counter()
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status)
        elapsed = time.perf_counter() - t
    if proc.returncode:
        raise RuntimeError(prefix)
    q, r = [
        json.loads((args.work / p / "preparation.json").read_text())
        for p in (query, receiver)
    ]
    rot, trans, fit, nca = preparation.alignment(q, r)
    refs = np.array(q["mapped"]) @ rot + trans
    key = dict(
        split="confirmation"
        if args.split == "confirmation"
        else ("development" if query.startswith("1FK") else "transfer"),
        target=q["target"],
        query=query,
        receiver=receiver,
        seed=seed,
        cap=cap,
        arm=arm,
    )
    poses = []
    for rank, block in enumerate(
        prefix.with_suffix(".pdbqt").read_text().split("ENDMDL")[:-1], 1
    ):
        lines = block.splitlines()
        xyz = {
            int(l[6:11]): preparation.coords(l)
            for l in lines
            if l.startswith(("ATOM", "HETATM"))
        }
        ordered = np.array([xyz[i] for i in q["serial"]])
        rmsd = np.sqrt(np.sum((ordered - refs) ** 2, axis=2).mean(axis=1)).min()
        score = float(
            next(l.split()[3] for l in lines if l.startswith("REMARK VINA RESULT:"))
        )
        poses.append(dict(**key, rank=rank, score=score, rmsd=rmsd))
    hits = [p["rmsd"] <= 2.0 for p in poses]
    counters = re.search(
        r"hunt_evals=(\d+).*?second_evals=(\d+)", prefix.with_suffix(".log").read_text()
    )
    return dict(
        **key,
        top1=int(hits[0]),
        top5=int(any(hits[:5])),
        oracle=int(any(hits)),
        top1_rmsd=poses[0]["rmsd"],
        pool_best=min(p["rmsd"] for p in poses),
        count=len(poses),
        wall_s=elapsed,
        cpu_s=usage.ru_utime + usage.ru_stime,
        peak_rss_KiB=usage.ru_maxrss,
        evals=sum(map(int, counters.groups())) if counters else "",
        pocket_ca=nca,
        ca_rmsd=fit,
    ), poses


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "out", "parent", "wall", "joint", "official"):
        p.add_argument(name, type=Path)
    p.add_argument(
        "--split",
        choices=("development", "transfer", "confirmation", "vina_reference"),
        required=True,
    )
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    a = p.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    (a.out / "runs").mkdir(parents=True, exist_ok=True)
    if a.split in ("development", "vina_reference"):
        pairs = list(itertools.product(("1FKG", "1FKH"), repeat=2))
    else:
        pairs = [
            (pdb, pdb)
            for _, pdb, _ in ENTRIES
            if json.loads((a.work / pdb / "preparation.json").read_text())["status"]
            == "ready"
        ]
    seeds, caps = (14591, 30553), (3000, 12000)
    arms = ("parent", "wall", "joint", "official")
    if a.split == "confirmation":
        seeds, caps, arms = (104729, 130363), (3000, 12000, 48000), (*arms, "vina")
    if a.split == "vina_reference":
        pairs += [
            (pdb, pdb)
            for _, pdb, _ in ENTRIES
            if json.loads((a.work / pdb / "preparation.json").read_text())["status"]
            == "ready"
        ]
        arms = ("vina",)
    jobs = [
        (*pair, seed, cap, arm)
        for pair, seed, cap, arm in itertools.product(pairs, seeds, caps, arms)
    ]
    random.Random(831337).shuffle(jobs)
    rows, poses = [], []
    for job in jobs:
        result, modes = run(a, job)
        rows.append(result)
        poses.extend(modes)
    for kind, data in (("cells", rows), ("models", poses)):
        write(a.out / f"{a.split}_{kind}.csv", data)
    hashes = {
        arm: hashlib.sha256(getattr(a, arm).read_bytes()).hexdigest()
        for arm in ("parent", "wall", "joint", "official")
    }
    (a.out / "binaries.json").write_text(json.dumps(hashes, indent=2) + "\n")


if __name__ == "__main__":
    main()
