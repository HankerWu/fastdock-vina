"""Verify measured cost of the reduced-budget closure arm on identical output pools."""

import argparse
import hashlib
import os
import random
import subprocess
import tempfile
import time
from pathlib import Path

from run import REFINE, TRAJECTORY, cases, read, write


def main(a):
    os.sched_setaffinity(0, {a.cpu})
    _, _, rec, lig, box = next(c for c in cases(a.work, "macro") if c[0] == "1FKB")
    reference = {
        (int(r["seed"]), r["arm"], int(r["cap"])): r for r in read(a.reference)
    }
    seeds = (458921, 740317, 191009)
    arms = [("parent", a.parent, 3000), ("candidate", a.candidate, 2000)]
    rows = []
    with tempfile.TemporaryDirectory() as temp:
        out = Path(temp) / "out.pdbqt"
        for block in range(8):
            order = list(seeds)
            random.Random(3011 + block).shuffle(order)
            for seed in order:
                for arm, binary, cap in arms if block % 2 == 0 else reversed(arms):
                    ref = reference[seed, arm, cap]
                    cmd = [
                        binary,
                        "--receptor",
                        rec,
                        "--ligand",
                        lig,
                        "--config",
                        box,
                        "--num_modes",
                        "5",
                        "--energy_range",
                        "100",
                        "--cpu",
                        "1",
                        "--exhaustiveness",
                        "1",
                        "--seed",
                        str(seed),
                        "--max_evals",
                        str(cap),
                        "--search_counters",
                        "--verbosity",
                        "0",
                        "--out",
                        out,
                    ]
                    start = time.perf_counter_ns()
                    p = subprocess.run(cmd, capture_output=True, text=True, check=True)
                    elapsed = time.perf_counter_ns() - start
                    log = p.stdout + p.stderr
                    sha = hashlib.sha256(out.read_bytes()).hexdigest()
                    fields = (
                        "hunt_calls",
                        "hunt_evals",
                        "second_calls",
                        "second_evals",
                        "explicit_calls",
                        "explicit_evals",
                    )
                    actual = (
                        *TRAJECTORY.search(log).groups(),
                        *REFINE.search(log).groups(),
                    )
                    if sha != ref["pool_sha256"] or actual != tuple(
                        ref[k] for k in fields
                    ):
                        raise ValueError((block, seed, arm))
                    rows.append(
                        dict(
                            block=block,
                            seed=seed,
                            arm=arm,
                            cap=cap,
                            elapsed_ns=elapsed,
                            joint=ref["joint"],
                            pool_sha256=sha,
                        )
                    )
    write(a.output, rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "parent", "candidate", "reference", "output"):
        p.add_argument(k, type=Path)
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    main(p.parse_args())
