"""Compare every upstream and monitored production output with clean executable replays."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

from backtracking_panel import CASES, paired


def replay(work, release, pristine):
    records = []
    m = {q: json.loads((work / q / "preparation.json").read_text()) for q in CASES}
    with (work / "cells.csv").open() as f:
        cells = [r for r in csv.DictReader(f) if r["arm"] in {"parent", "upstream"}]
    for r in cells:
        arm = r["arm"]
        label, binary = (
            ("release", release) if arm == "upstream" else ("pristine", pristine)
        )
        other, _ = paired.execute(
            work,
            m,
            r["query"],
            r["receiver"],
            int(r["seed"]),
            int(r["cap"]),
            label,
            binary,
            False,
        )
        ref = (work / "runs" / (r["key"] + ".pdbqt")).read_bytes()
        actual = (work / "runs" / (other["key"] + ".pdbqt")).read_bytes()
        if actual != ref:
            raise ValueError(r["key"] + " coordinate parity failure")
        if arm == "parent":
            assert all(
                int(r[k]) == other[k]
                for k in (
                    "hunt_calls",
                    "hunt_evals",
                    "second_calls",
                    "second_evals",
                    "explicit_calls",
                    "explicit_evals",
                )
            )
        records.append(
            {
                "key": r["key"],
                "replay": label,
                "pdbqt_sha256": hashlib.sha256(ref).hexdigest(),
                "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
            }
        )
    paired.write(work / "parity.csv", records)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument("--release", type=Path, required=True)
    ap.add_argument("--pristine", type=Path, required=True)
    ap.add_argument("--cpu", type=int, required=True)
    args = ap.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    replay(args.work, args.release, args.pristine)
