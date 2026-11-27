"""Repeat discovery from archived structures; write paired native outputs and provenance."""

import argparse
import hashlib
import itertools
import json
import os
import random
import re
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parents[1] / "pair_gradient"
sys.path.insert(0, str(PRIOR / "pose"))
import paired


def read_files(archive, work, ids):
    with tarfile.open(archive) as tar:
        files = [
            m
            for m in tar
            if len(Path(m.name).parts) == 2
            and Path(m.name).parts[0] in ids
            and Path(m.name).name
            in ("receptor.pdbqt", "ligand.pdbqt", "preparation.json")
        ]
        tar.extractall(work, members=files, filter="data")


def main(args):
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    (work / "runs").mkdir(exist_ok=True)
    # The paired scorer uses original deposited REN sequences to fix author numbering.
    with tarfile.open(PRIOR / "pose/inputs.tar.xz") as archive:
        raw = [
            m for m in archive if m.name.startswith("raw/") and m.name.endswith(".pdb")
        ]
        archive.extractall(work, members=raw, filter="data")
    families = (("3EML", "3UZA"), ("3G6Z", "3SFC"), ("1E3G", "2AM9"))
    read_files(PRIOR / "pose/inputs.tar.xz", work, set(itertools.chain(*families[:2])))
    read_files(PRIOR / "holdout/inputs.tar.xz", work, set(families[-1]))
    meta = {
        p: json.loads((work / p / "preparation.json").read_text())
        for p in itertools.chain(*families)
    }
    paired.renin_correspondence(work, meta)
    arms = ("parent", "candidate", "ungated")
    todo = [
        (s, r, seed, cap)
        for family in families
        for s, r in itertools.product(family, repeat=2)
        for seed, cap in itertools.product((75467, 98251), (3000, 12000))
    ]
    rng = random.Random(270113)
    rng.shuffle(todo)
    os.sched_setaffinity(
        0, {args.cpu if args.cpu is not None else min(os.sched_getaffinity(0))}
    )
    cells, modes = [], []
    for s, r, seed, cap in todo:
        order = list(arms)
        rng.shuffle(order)
        case = {}
        for arm in order:
            res = paired.execute(
                work, meta, s, r, seed, cap, arm, getattr(args, arm), False
            )
            row, poses = res
            key = row["key"]
            log = (work / "runs" / (key + ".log")).read_text()
            stats = dict(
                re.findall(
                    r"(\w+)=([\w.\-+]+)",
                    next((t for t in log.splitlines() if "Deferred hunt:" in t), ""),
                )
            )
            if arm != "parent" and len(stats) != 9:
                raise ValueError(f"Missing retention counters {key} {stats}")
            for col in (
                "skipped",
                "checks",
                "records",
                "eligible",
                "replaced",
                "slots",
                "slot",
                "pending_energy",
                "previous_energy",
            ):
                row[col] = stats.get(col, "")
            row["pool_sha256"] = hashlib.sha256(
                (work / "runs" / (key + ".pdbqt")).read_bytes()
            ).hexdigest()
            cells.append(row)
            modes.extend(poses)
            case[arm] = row
        invariant = (
            "hunt_calls",
            "hunt_evals",
            "second_calls",
            "second_evals",
            "explicit_inputs",
        )
        for arm in arms[1:]:
            for name in invariant:
                if case[arm][name] != case["parent"][name]:
                    raise ValueError((s, r, seed, arm, name))
            if (
                case[arm]["replaced"] == "0"
                and case[arm]["pool_sha256"] != case["parent"]["pool_sha256"]
            ):
                raise ValueError((s, r, seed, arm, "unchanged pool differs"))
            if case[arm]["slots"] != str(case["parent"]["explicit_inputs"]):
                raise ValueError("slot count")
        if (
            case["candidate"]["replaced"] == "1"
            and case["candidate"]["pool_sha256"] != case["ungated"]["pool_sha256"]
        ):
            raise ValueError("two identical interventions differ")
    paired.write(work / "cells.csv", cells)
    paired.write(work / "modes.csv", modes)
    (work / "builds.json").write_text(
        json.dumps(
            {
                a: hashlib.sha256(getattr(args, a).read_bytes()).hexdigest()
                for a in arms
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    for arm in ("parent", "candidate", "ungated"):
        p.add_argument("--" + arm, type=Path, required=True)
    p.add_argument("--cpu", type=int)
    main(p.parse_args())
