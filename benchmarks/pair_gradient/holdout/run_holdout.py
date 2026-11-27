"""Five-arm paired docking and optional passive, byte-checked stage/cost replays."""

import argparse
import hashlib
import itertools
import json
import random
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pose"))
import paired

ARMS = ("upstream", "parent", "linear", "hermite", "shift")
TRACE = re.compile(r"^\**STAGE (\w+) (\d+)\n(.*?)^END STAGE$", re.M | re.S)


def execute(binary, common, path, counters):
    cmd = [str(binary.resolve()), *common, "--out", str(path)]
    if counters:
        cmd.append("--search_counters")
    start = time.perf_counter()
    p = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return p.stdout + p.stderr, time.perf_counter() - start


def stage_data(trace, key, meta, receiver):
    rot, trans = paired.base.alignment(meta, receiver)[:2]
    ref = np.asarray(meta["mapped"]) @ rot + trans
    rows = []
    for stage, idx, block in TRACE.findall(trace):
        exact, pdbqt = block.split("PDBQT\n", 1)
        atoms = paired.atoms(pdbqt.splitlines())
        xyz = np.array([atoms[i][0] for i in meta["serial"]])
        rms = np.sqrt(((ref - xyz) ** 2).sum(axis=2).mean(axis=1)).min()
        rows.append(
            dict(
                key=key,
                stage=stage,
                input_index=int(idx),
                rmsd_A=float(rms),
                energy=float.fromhex(exact.splitlines()[0]),
                coords_sha256=hashlib.sha256(
                    exact.partition("\n")[2].encode()
                ).hexdigest(),
                sha256=hashlib.sha256(exact.encode()).hexdigest(),
            )
        )
    if not rows:
        raise ValueError("Missing stage trace: " + key)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    for arm in ARMS:
        p.add_argument("--" + arm, type=Path, required=True)
    p.add_argument("--counts-replay", type=Path, required=True)
    p.add_argument("--trace-linear", type=Path)
    p.add_argument("--trace-shift", type=Path)
    p.add_argument(
        "--replay-only",
        action="store_true",
        help="reuse production timings and PDBQT outputs",
    )
    args = p.parse_args()
    work = args.work.resolve()
    meta = {
        code: json.loads((work / code / "preparation.json").read_text())
        for code in json.loads((work / "cases.json").read_text())
    }
    specs = list(
        itertools.product(
            meta,
            [r for r in meta if meta[r]["receptor_status"] == "ready"],
            paired.SEEDS,
            paired.CAPS,
        )
    )
    random.Random(20261127).shuffle(specs)
    directory = work / "runs"
    directory.mkdir(exist_ok=True)
    cells, modes, stages = [], [], []
    for i, (q, r, seed, cap) in enumerate(specs):
        start_arm = i % len(ARMS)
        for arm in ARMS[start_arm:] + ARMS[:start_arm]:
            key = f"{q}_{r}_{seed}_{cap}_{arm}"
            out = directory / (key + ".pdbqt")
            logpath = directory / (key + ".log")
            timepath = directory / (key + ".time")
            common = [
                "--receptor",
                str(work / r / "receptor.pdbqt"),
                "--ligand",
                str(work / q / "ligand.pdbqt"),
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
            for axis, v in zip("xyz", meta[r]["center"]):
                common += ["--center_" + axis, str(v), "--size_" + axis, "24"]
            if not args.replay_only:
                log, wall = execute(getattr(args, arm), common, out, arm != "upstream")
                logpath.write_text(log)
                timepath.write_text(str(wall) + "\n")
            log = logpath.read_text()
            diag = (
                args.counts_replay
                if arm == "upstream"
                else getattr(args, "trace_" + arm, None)
            )
            parity = ""
            if diag:
                sidecar = directory / (
                    key + (".counts" if arm == "upstream" else ".stages")
                )
                if not sidecar.exists() or not args.replay_only:
                    replay = directory / (key + ".replay")
                    diaglog, _ = execute(diag, common, replay, arm != "upstream")
                    if replay.read_bytes() != out.read_bytes():
                        raise RuntimeError(
                            "Instrumentation changed pose output: " + key
                        )
                    sidecar.write_text(diaglog)
                    replay.unlink()
                parity = "identical"
                log = sidecar.read_text()
                if arm != "upstream":
                    stages.extend(stage_data(log, key, meta[q], meta[r]))
            row = dict(
                key=key,
                target="AR",
                query=q,
                receiver=r,
                seed=seed,
                cap=cap,
                arm=arm,
                kind="cognate" if q == r else "cross",
                seconds=float(timepath.read_text()),
                diagnostic_replay=parity,
                **paired.counts(
                    log.replace("Search trajectory:", "\nSearch trajectory:")
                ),
            )
            if any(
                row[k] == "" for k in ("hunt_evals", "explicit_evals", "second_evals")
            ):
                raise ValueError("Missing joint evaluation data: " + key)
            modes.extend(paired.inspect(directory, meta[q], meta[r], row))
            cells.append(row)
        print(q, r, seed, cap, flush=True)
    paired.write(work / "cells.csv", cells)
    paired.write(work / "modes.csv", modes)
    if stages:
        paired.write(work / "stages.csv", stages)


if __name__ == "__main__":
    main()
