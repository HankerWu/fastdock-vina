"""Run the two frozen historical structures. Score printed modes by two atom maps."""

import argparse
import hashlib
import itertools
import json
import os
import random
import re
import subprocess
import sys
import time
from pathlib import Path

EXTERNAL = Path(__file__).resolve().parents[1] / "macrocycle_proposal" / "external"
sys.path.insert(0, str(EXTERNAL))
from external_run import write
from external_score import main as score


def main(args):
    os.sched_setaffinity(0, {args.cpu})
    raw = args.output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    cases = list(
        itertools.product(("5T9W", "7MBO"), (458921, 740317, 191009), (3000, 12000))
    )
    random.Random(36).shuffle(cases)
    records = []
    arms = ("upstream", "parent", "always", "adaptive")
    for j, (pdb, seed, cap) in enumerate(cases):
        directory = args.inputs / pdb
        meta = json.loads((directory / "preparation.json").read_text())
        box = [
            v
            for d, c in zip("xyz", meta["center"])
            for v in ("--center_" + d, str(c), "--size_" + d, "26")
        ]
        for arm in arms[j % 4 :] + arms[: j % 4]:
            key = f"{pdb}_{seed}_{cap}_{arm}"
            dest = raw / (key + ".pdbqt")
            cmd = [
                "--ligand",
                str(directory / "ligand.pdbqt"),
                "--receptor",
                str(directory / "receptor.pdbqt"),
                *box,
                "--cpu",
                "1",
                "--seed",
                str(seed),
                "--exhaustiveness",
                "1",
                "--num_modes",
                "5",
                "--max_evals",
                str(cap),
                "--energy_range",
                "100",
                "--verbosity",
                "0",
                "--out",
                str(dest),
            ]
            opts = [] if arm == "upstream" else ["--search_counters"]
            start = time.perf_counter_ns()
            with (raw / (key + ".log")).open("w") as log:
                proc = subprocess.Popen(
                    [getattr(args, arm), *cmd, *opts],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                )
                _, status, usage = os.wait4(proc.pid, 0)
                proc.returncode = os.waitstatus_to_exitcode(status)
            elapsed = time.perf_counter_ns() - start
            if proc.returncode:
                raise RuntimeError(key)
            log = (raw / (key + ".log")).read_text()
            if arm == "upstream":
                replay = raw / (key + ".counter-pose")
                cmd[-1] = str(replay)
                log = subprocess.check_output(
                    [args.upstream_counts, *cmd], stderr=subprocess.STDOUT, text=True
                )
                (raw / (key + ".counts")).write_text(log)
                if replay.read_bytes() != dest.read_bytes():
                    raise ValueError("Counter instrumentation changed " + key)
                replay.unlink()
            counts = {
                k: int(v)
                for k, v in re.findall(
                    r"(hunt_calls|hunt_evals|second_calls|second_evals|explicit_calls|explicit_evals)=(\d+)",
                    log,
                )
            }
            if len(counts) != 6:
                raise ValueError("Missing measured counts for " + key)
            flips = re.search(r"rigid_reorientations=(\d+)", log)
            records.append(
                dict(
                    key=key,
                    pdb=pdb,
                    seed=seed,
                    cap=cap,
                    arm=arm,
                    elapsed_ns=elapsed,
                    cpu_s=usage.ru_utime + usage.ru_stime,
                    peak_rss_kib=usage.ru_maxrss,
                    flips=int(flips[1]) if flips else 0,
                    **counts,
                    joint=sum(v for k, v in counts.items() if k.endswith("evals")),
                    modes=dest.read_text().count("ENDMDL"),
                    sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),
                )
            )
    write(args.output / "runs.csv", records)
    score(args.inputs, args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for field in (
        "inputs",
        "output",
        "parent",
        "upstream",
        "upstream-counts",
        "always",
        "adaptive",
    ):
        parser.add_argument("--" + field, type=Path, required=True)
    parser.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    main(parser.parse_args())
