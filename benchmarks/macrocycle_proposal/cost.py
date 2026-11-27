"""Counterbalance identical production predictions across four independent process launches."""

import argparse
import itertools
import os
from pathlib import Path
import random
import resource
import subprocess
import tempfile
import time

from panel import cases, command, counts, read, sha, write


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "results", "out", "parent", "reset", "ccd"):
        p.add_argument(name, type=Path)
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    a = p.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    reference = {r["key"]: r for r in read(a.results / "runs.csv")}
    panel = list(cases(a.inputs, "macro"))
    cells = list(
        itertools.product(panel, (740317, 458921), ((3000, 2000), (12000, 9000)))
    )
    rows = []
    with tempfile.TemporaryDirectory() as temp:
        pose = Path(temp) / "pose.pdbqt"
        for block in range(4):
            order = cells.copy()
            random.Random(101 + block // 2).shuffle(order)
            for i, (case, seed, (cap, lower)) in enumerate(order):
                arms = [
                    ("parent", cap),
                    ("reset", cap),
                    ("ccd", cap),
                    ("reset", lower),
                    ("ccd", lower),
                ]
                j = (i + block // 2) % len(arms)
                arms = arms[j:] + arms[:j]
                if block % 2:
                    arms.reverse()
                for arm, actual in arms:
                    key = f"{case[0]}_{case[1]}_{seed}_{actual}_{arm}"
                    r = reference[key]
                    cmd = command(case, seed, actual, getattr(a, arm), pose)
                    start_cpu = resource.getrusage(resource.RUSAGE_CHILDREN)
                    start = time.perf_counter_ns()
                    log = subprocess.check_output(
                        cmd, text=True, stderr=subprocess.STDOUT
                    )
                    elapsed = time.perf_counter_ns() - start
                    end_cpu = resource.getrusage(resource.RUSAGE_CHILDREN)
                    cpu_seconds = (end_cpu.ru_utime + end_cpu.ru_stime) - (
                        start_cpu.ru_utime + start_cpu.ru_stime
                    )
                    if r["pool_sha256"] != sha(pose) or any(
                        int(r[k]) != v for k, v in counts(log).items()
                    ):
                        raise ValueError(key)
                    rows.append(
                        dict(
                            block=block,
                            key=key,
                            joint=r["joint"],
                            elapsed_ns=elapsed,
                            cpu_seconds=round(cpu_seconds, 6),
                        )
                    )
    write(a.out, rows)


if __name__ == "__main__":
    main()
