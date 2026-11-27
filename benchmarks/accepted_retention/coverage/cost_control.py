"""Cost sensitivity of all cells where the new selector changes a 2 Å endpoint."""

import argparse
import csv
import hashlib
import os
import random
import re
import subprocess
import time
from pathlib import Path

from probe import command, counters, mode_distances, references, table


def main(a):
    os.sched_setaffinity(0, {a.cpu})
    rows = list(csv.DictReader((a.results / "blind.csv").open()))
    indexed = {(r["key"], r["arm"]): r for r in rows}
    changes = [
        r
        for r in rows
        if r["arm"] == "parent"
        and any(
            (float(r[k]) <= 2) != (float(indexed[r["key"], "gated"][k]) <= 2)
            for k in ("top1", "top5", "oracle")
        )
    ]
    measurements = []
    used = {}
    rng = random.Random(701723)
    for r in changes:
        serial, refs = references(a.work, r)
        for repeat in range(7):
            arms = [(arm, cap) for arm in ("parent", "gated") for cap in (11900, 12000)]
            rng.shuffle(arms)
            for arm, cap in arms:
                out = a.work / "cost_replay.pdbqt"
                cmd = command(a.work, r, a.parent if arm == "parent" else a.gated, out)
                cmd[cmd.index("--max_evals") + 1] = str(cap)
                start = time.perf_counter_ns()
                p = subprocess.run(cmd, capture_output=True, check=True)
                ns = time.perf_counter_ns() - start
                sha = hashlib.sha256(out.read_bytes()).hexdigest()
                key = r["key"], arm, cap
                costs = dict(re.findall(rb"(\w+)=(\d+)", b" ".join(counters(p.stdout))))
                joint = sum(
                    int(costs[x])
                    for x in (b"hunt_evals", b"second_evals", b"explicit_evals")
                )
                if used.setdefault(key, (sha, joint)) != (sha, joint) or (
                    cap == 12000 and sha != indexed[r["key"], arm]["sha256"]
                ):
                    raise ValueError("Deterministic replay failed")
                score = mode_distances(out.read_bytes(), serial, refs)
                measurements.append(
                    {
                        "key": r["key"],
                        "arm": arm,
                        "cap": cap,
                        "repeat": repeat,
                        "ns": ns,
                        "joint": joint,
                        "top1": score[0],
                        "top5": min(score[:5]),
                        "oracle": min(score),
                        "sha256": sha,
                    }
                )
    if measurements:
        table(a.results / "cost_control.csv", measurements)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "results", "parent", "gated"):
        p.add_argument(k, type=Path)
    p.add_argument("--cpu", type=int, required=True)
    main(p.parse_args())
