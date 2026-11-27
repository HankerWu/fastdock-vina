"""Run the historical macrocycles with identical inputs and one trajectory per prediction."""

import argparse
import hashlib
import itertools
import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH / "bond_distance"))
from run import cases, write, read, TRAJECTORY, REFINE
from analyze import structural


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command(case, seed, cap, binary, pose, counters=True):
    pdb, cid, rec, lig, box = case
    args = [
        binary,
        "--receptor",
        rec,
        "--ligand",
        lig,
        "--config",
        box,
        "--cpu",
        "1",
        "--exhaustiveness",
        "1",
        "--seed",
        str(seed),
        "--max_evals",
        str(cap),
        "--num_modes",
        "5",
        "--energy_range",
        "100",
        "--verbosity",
        "0",
        "--out",
        pose,
    ]
    return args + (["--search_counters"] if counters else [])


def counts(log):
    t, e = TRAJECTORY.search(log), REFINE.search(log)
    if not t or not e:
        raise ValueError("Missing trajectory or refinement counters")
    return dict(
        zip(
            (
                "hunt_calls",
                "hunt_evals",
                "second_calls",
                "second_evals",
                "explicit_calls",
                "explicit_evals",
            ),
            map(int, t.groups() + e.groups()),
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "inputs",
        "out",
        "upstream",
        "upstream_counts",
        "parent",
        "reset",
        "ccd",
    ):
        parser.add_argument(name, type=Path)
    parser.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    a = parser.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    raw = a.out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    panel = list(cases(a.inputs, "macro"))
    cells = list(itertools.product(panel, (740317, 458921), (2000, 3000, 9000, 12000)))
    random.Random(458921).shuffle(cells)
    rows = []
    for i, (case, seed, cap) in enumerate(cells):
        arms = (
            ["upstream", "parent", "ccd", "reset"]
            if cap in (3000, 12000)
            else ["ccd", "reset"]
        )
        for arm in arms[i % len(arms) :] + arms[: i % len(arms)]:
            pdb, cid, rec, lig, box = case
            key = f"{pdb}_{cid}_{seed}_{cap}_{arm}"
            pose = raw / (key + ".pdbqt")
            args = command(case, seed, cap, getattr(a, arm), pose, arm != "upstream")
            start = time.perf_counter_ns()
            log = subprocess.check_output(args, text=True, stderr=subprocess.STDOUT)
            elapsed = time.perf_counter_ns() - start
            (raw / (key + ".log")).write_text(log)
            if arm == "upstream":
                replay = raw / (key + ".replay")
                log = subprocess.check_output(
                    command(case, seed, cap, a.upstream_counts, replay, False),
                    text=True,
                    stderr=subprocess.STDOUT,
                )
                if replay.read_bytes() != pose.read_bytes():
                    raise ValueError(key + ": instrumentation changed output")
                replay.unlink()
                (raw / (key + ".counts")).write_text(log)
            metrics = counts(log)
            scores = [
                s.split()[3]
                for s in pose.read_text().splitlines()
                if s.startswith("REMARK VINA RESULT:")
            ]
            rows.append(
                dict(
                    key=key,
                    pdb=pdb,
                    compound=cid,
                    seed=seed,
                    cap=cap,
                    arm=arm,
                    elapsed_ns=elapsed,
                    first_score=scores[0],
                    modes=len(scores),
                    **metrics,
                    joint=sum(v for k, v in metrics.items() if k.endswith("_evals")),
                    pool_sha256=sha(pose),
                )
            )
    write(a.out / "runs.csv", rows)
    structural(a.inputs, a.out, read(a.out / "runs.csv"))
    (a.out / "binaries.json").write_text(
        json.dumps(
            {
                k: sha(getattr(a, k))
                for k in ("upstream", "upstream_counts", "parent", "reset", "ccd")
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
