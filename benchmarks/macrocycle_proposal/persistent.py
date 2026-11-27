"""Pair the persistent constraint with fresh parent processes on the historical panel."""

import argparse
import itertools
import os
from pathlib import Path
import random
import subprocess
import time
from panel import cases, command, counts, read, sha, structural, write


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "previous", "out", "parent", "lock"):
        p.add_argument(name, type=Path)
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    a = p.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    panel = list(cases(a.inputs, "macro"))
    historical = {r["key"]: r for r in read(a.previous / "runs.csv")}
    raw = a.out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    cells = list(
        itertools.product(panel, (740317, 458921), ((3000, 2000), (12000, 9000)))
    )
    rows = []
    times = []
    for block in range(4):
        order = cells.copy()
        random.Random(511 + block // 2).shuffle(order)
        for i, (case, seed, (cap, low)) in enumerate(order):
            pdb, cid = case[:2]
            arms = [("parent", cap), ("lock", cap), ("lock", low)]
            j = i % 3
            arms = arms[j:] + arms[:j]
            if block % 2:
                arms.reverse()
            for arm, budget in arms:
                key = f"{pdb}_{cid}_{seed}_{budget}_{arm}"
                pose = raw / (key + ".pdbqt")
                cmd = command(case, seed, budget, getattr(a, arm), pose)
                start = time.perf_counter_ns()
                log = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
                elapsed = time.perf_counter_ns() - start
                (raw / (key + ".log")).write_text(log)
                scores = [
                    s.split()[3]
                    for s in pose.read_text().splitlines()
                    if s.startswith("REMARK VINA RESULT:")
                ]
                counter = counts(log)
                digest = sha(pose)
                joint = sum(v for k, v in counter.items() if k.endswith("_evals"))
                item = dict(
                    key=key,
                    pdb=pdb,
                    compound=cid,
                    seed=seed,
                    cap=budget,
                    arm=arm,
                    elapsed_ns=elapsed,
                    first_score=scores[0],
                    modes=len(scores),
                    **counter,
                    joint=joint,
                    pool_sha256=digest,
                )
                if block == 0:
                    rows.append(item)
                    if arm == "lock":
                        historical[key] = item
                prior = historical[key]
                if digest != prior["pool_sha256"] or any(
                    v != int(prior[k]) for k, v in counter.items()
                ):
                    raise ValueError("Non-reproducible prediction: " + key)
                times.append(
                    dict(block=block, key=key, elapsed_ns=elapsed, joint=joint)
                )
    write(a.out / "runs.csv", rows)
    write(a.out / "timing.csv", times)
    structural(a.inputs, a.out, rows)


if __name__ == "__main__":
    main()
