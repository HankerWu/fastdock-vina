"""Recreate shared inputs or compare registered binaries on all eligible structures."""

import argparse
import hashlib
import itertools
import json
import os
import random
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
sys.path.insert(0, str(BENCH / "rigid_flip"))
import holdout

CASES = {
    "2LBD": ("RARG", "REA", "A"),
    "1FBY": ("RXRA", "9CR", "A"),
    "1Q41": ("GSK3B", "IXM", "A"),
    "1Q5K": ("GSK3B", "TMU", "A"),
}
SEEDS = (51871, 67213)
CAPS = (3000, 12000, 2550, 10200)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "analyze"))
    parser.add_argument("work", type=Path)
    for name in ("parent", "upstream", "counts", "always", "adaptive"):
        parser.add_argument("--" + name, type=Path)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    if args.action == "prepare":
        sys.path.insert(0, str(BENCH / "archive_overlap/transfer"))
        import prepare_inputs

        prepare_inputs.main(work, CASES)
        return
    meta = {p: json.loads((work / p / "preparation.json").read_text()) for p in CASES}
    arms = ("upstream", "parent", "always", "adaptive")
    directory = work / "runs"
    directory.mkdir(exist_ok=True)
    cells, excluded = [], []
    for q, r in itertools.product(sorted(CASES), repeat=2):
        if CASES[q][0] != CASES[r][0]:
            continue
        if meta[r]["receptor_status"] != "ready":
            excluded.append([q, r, meta[r]["receptor_status"]])
            continue
        try:
            holdout.paired.base.alignment(meta[q], meta[r])
        except ValueError as exc:
            excluded.append([q, r, str(exc)])
            continue
        cells.extend((q, r, s, c) for s, c in itertools.product(SEEDS, CAPS))
    (work / "exclusions.json").write_text(json.dumps(excluded, indent=2) + "\n")
    if args.action == "run":
        binaries = {arm: getattr(args, arm) for arm in (*arms, "counts")}
        if not all(binaries.values()):
            parser.error("run requires all four arms and the upstream counter replay")
        builds = {
            arm: hashlib.sha256(p.read_bytes()).hexdigest()
            for arm, p in binaries.items()
        }
        build_record = work / "builds.json"
        if build_record.exists() and json.loads(build_record.read_text()) != builds:
            raise ValueError(
                "existing runs belong to different binaries; use a fresh work directory"
            )
        build_record.write_text(json.dumps(builds, indent=2) + "\n")
        order = random.Random(112711)
        order.shuffle(cells)
        holdout.CASES = CASES
        for cell in cells:
            permuted = list(
                arms if cell[3] in (3000, 12000) else ("always", "adaptive")
            )
            order.shuffle(permuted)
            for arm in permuted:
                holdout.run(work, binaries[arm], args.counts, meta, *cell, arm)
        return
    rows, modes = [], []
    for q, r, s, c in sorted(cells):
        for arm in arms if c in (3000, 12000) else ("always", "adaptive"):
            key = f"{q}_{r}_{s}_{c}_{arm}"
            row = json.loads((directory / (key + ".json")).read_text())
            del row["rigid_reorientations"]
            log = (directory / (key + ".log")).read_text()
            for k in ("axes", "other", "offered", "redirected"):
                row["silent_" + k] = ""
            for line in log.splitlines():
                if "Silent torsions:" in line:
                    for k, v in re.findall(r"(\w+)=(\d+)", line):
                        row["silent_" + k] = int(v)
            modes.extend(holdout.paired.inspect(directory, meta[q], meta[r], row))
            rows.append(row)
    holdout.paired.write(work / "cells.csv", rows)
    holdout.paired.write(work / "modes.csv", modes)


if __name__ == "__main__":
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    main()
