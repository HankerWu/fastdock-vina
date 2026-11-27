"""Run every four-trajectory Gal3 pair and compare each PDBQT to a rank-stable offline filter."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import random
import re
import subprocess
import sys
import tarfile

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]


def models(data):
    result = []
    for part in data.decode().split("ENDMDL")[:-1]:
        lines = part.strip().splitlines()
        if not lines[0].startswith("MODEL"):
            raise ValueError("Missing model boundary")
        body = "\n".join(lines[1:])
        coords = np.array(
            [
                [float(l[i : i + 8]) for i in (30, 38, 46)]
                for l in lines
                if l.startswith(("ATOM  ", "HETATM"))
                and l.split()[-1] not in ("H", "HS", "HD")
            ]
        )
        result.append((body, coords))
    return result


def selection(pool):
    chosen = []
    nearest = 1e30
    ambiguous = 0
    for i, (_, xyz) in enumerate(pool):
        ds = [np.sqrt(np.square(xyz - pool[j][1]).sum(axis=1).mean()) for j in chosen]
        if ds:
            nearest = min(nearest, min(ds))
            ambiguous += sum(abs(d - 1) <= np.sqrt(3) * 0.001 for d in ds)
        if all(d >= 1 for d in ds):
            chosen.append(i)
    return chosen, nearest, ambiguous


def archived(name):
    with tarfile.open(name) as t:
        for entry in t:
            if entry.name.endswith(".pdbqt"):
                yield entry.name, t.extractfile(entry).read()


def offline():
    panels = {
        "development": ("development/raw.tar.xz", "_fixed.pdbqt"),
        "affinity_one": ("affinity/single.tar.xz", "_candidate_"),
        "affinity_four": ("affinity/four.tar.xz", "_candidate_"),
        "affinity_reduced": ("affinity/reduced.tar.xz", "_candidate_"),
        "transfer": ("transfer/outputs.tar.xz", "_candidate.pdbqt"),
        "transfer_challenge": ("transfer/refill/outputs.tar.xz", "_candidate.pdbqt"),
        "transfer_reduced": ("transfer/mechanism/runs.tar.xz", "_candidate.pdbqt"),
    }
    for panel, (name, suffix) in panels.items():
        for path, data in archived(HERE.parent / name):
            if suffix not in path or (
                panel == "transfer_reduced" and not path.startswith("reduced/")
            ):
                continue
            pool = models(data)
            selected, mind, amb = selection(pool)
            yield dict(
                panel=panel,
                condition=Path(path).stem,
                modes=len(pool),
                kept=len(selected),
                min_to_retained_A="" if len(pool) < 2 else mind,
                rounding_ambiguous=amb,
                removed=",".join(
                    str(i + 1) for i in range(len(pool)) if i not in selected
                ),
            )


def output(path, rows):
    if not rows:
        raise ValueError("Empty output")
    with path.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(args):
    args.out.mkdir(parents=True, exist_ok=True)
    output(args.out / "offline.csv", list(offline()))
    cases = {
        Path(f).stem: x
        for f, x in archived(HERE.parent / "affinity/four.tar.xz")
        if "_candidate_" in f
    }
    keys = sorted(cases)
    random.Random(8753).shuffle(keys)
    rows = []
    raw = args.out / "raw"
    raw.mkdir(exist_ok=True)
    os.sched_setaffinity(0, {args.cpu})
    for i, key in enumerate(keys):
        _, r, q, seed, _, cap = key.split("_")
        expected = models(cases[key])
        indexes, _, amb = selection(expected)
        if amb:
            raise ValueError(key + " needs full-precision control")
        for arm, binary in [("old", args.old), ("new", args.new)][
            :: (-1 if i % 2 else 1)
        ]:
            name = raw / (key + "_" + arm)
            pose = name.with_suffix(".pdbqt")
            log = name.with_suffix(".log")
            cost = name.with_suffix(".json")
            args_cli = [
                "--config",
                args.inputs / (r + ".box"),
                "--ligand",
                args.inputs / (q + ".pdbqt"),
                "--receptor",
                args.inputs / (r + ".pdbqt"),
                "--cpu",
                "1",
                "--exhaustiveness",
                "4",
                "--max_evals",
                cap,
                "--seed",
                seed,
                "--num_modes",
                "9",
                "--energy_range",
                "100",
                "--verbosity",
                "1",
                "--search_counters",
                "--out",
                pose,
            ]
            subprocess.run(
                [
                    sys.executable,
                    "-S",
                    BENCH / "rigid_flip/process_cost.py",
                    log,
                    cost,
                    binary,
                    *args_cli,
                ],
                check=True,
            )
            actual = models(pose.read_bytes())
            ref = expected if arm == "old" else [expected[j] for j in indexes]
            if [p[0] for p in actual] != [p[0] for p in ref]:
                raise ValueError(key + " mode body mismatch")
            table = re.findall(
                r"^\s*(\d+)\s+([-+\d.eE]+)\s+[-+\d.eE]+\s+[-+\d.eE]+\s*$",
                log.read_text(),
                re.M,
            )
            if [int(rank) for rank, _ in table] != list(range(1, len(table) + 1)):
                raise ValueError(key + " non-contiguous console ranks")
            # The convergence table covers all energies; the writer applies --energy_range 100.
            eligible = [r for r, e in table if float(e) <= float(table[0][1]) + 100]
            if len(eligible) != len(ref):
                raise ValueError(key + " console/output mismatch")
            counts = "\n".join(
                x
                for x in log.read_text().splitlines()
                if x.startswith(
                    ("Search trajectory:", "Search counters:", "Hunt chain:")
                )
            )
            rows.append(
                dict(
                    condition=key,
                    arm=arm,
                    reported=len(actual),
                    sha256=hashlib.sha256(pose.read_bytes()).hexdigest(),
                    counters_sha256=hashlib.sha256(counts.encode()).hexdigest(),
                    **json.loads(cost.read_text()),
                )
            )
        if rows[-1]["counters_sha256"] != rows[-2]["counters_sha256"]:
            raise ValueError("changed search " + key)
    output(args.out / "cells.csv", rows)
    (args.out / "binaries.json").write_text(
        json.dumps(
            {
                k: hashlib.sha256(v.read_bytes()).hexdigest()
                for k, v in [("old", args.old), ("new", args.new)]
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("inputs", "out", "old", "new"):
        p.add_argument(k, type=Path)
    p.add_argument("--cpu", type=int, required=True)
    run(p.parse_args())
