"""Paired, counterbalanced end-to-end docking of the frozen Kd series."""

import argparse
import hashlib
import json
import os
import platform
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "accepted_retention/affinity"))
from run import COUNT, dock, read, write

SEED = 458921


def main(args):
    work = args.work.resolve()
    output = work / ("cost_check" if args.cost_check else "experiments")
    seed = 703109 if args.cost_check else SEED
    caps = (10800, 12000) if args.cost_check else (3000, 12000)
    raw = output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    os.sched_setaffinity(0, {args.cpu})
    pre = work / "prepared"
    arms = {
        name: getattr(args, name).resolve()
        for name in ("upstream", "parent", "halving", "quadratic")
    }
    if args.cost_check:
        del arms["upstream"]
    meta = {
        "seed": seed,
        "caps": caps,
        "arms": list(arms),
        "cpu": args.cpu,
        "platform": platform.platform(),
        "python": sys.version,
        "executables": {
            k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in arms.items()
        },
    }
    cells = []
    for row in read(work / "panel.csv"):
        rec, cid = row["pdb"], row["cid"]
        cell = (
            row["target"],
            rec,
            cid,
            pre / (rec + ".pdbqt"),
            pre / (rec + ".box"),
            pre / (cid + ".pdbqt"),
        )
        cells.extend((cell, cap) for cap in caps)
    random.Random(seed).shuffle(cells)
    result = []
    for i, (cell, cap) in enumerate(cells):
        order = list(arms)
        j = i % len(order)
        for arm in order[j:] + order[:j]:
            row = dock(cell, output, arms[arm], arm, seed, cap, arm != "upstream")
            if arm == "upstream":
                counted = dock(
                    cell, output, args.upstream_counts, arm, seed, cap, False, diag=True
                )
                f = raw / row["key"]
                if (
                    f.with_suffix(".pdbqt").read_bytes()
                    != f.with_suffix(".replay").read_bytes()
                ):
                    raise ValueError(
                        "Instrumented upstream output differs: " + row["key"]
                    )
                f.with_suffix(".replay").unlink()
                row.update({k: counted[k] for k in (*COUNT, "joint")})
            if not row["joint"]:
                raise ValueError("Missing counters: " + row["key"])
            result.append(row)
    write(output / "runs.csv", result)
    (output / "setup.json").write_text(json.dumps(meta, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    for opt in ("upstream", "upstream-counts", "parent", "halving", "quadratic"):
        parser.add_argument("--" + opt, required=True, type=Path)
    parser.add_argument("--cost-check", action="store_true")
    parser.add_argument("--cpu", type=int, required=True)
    main(parser.parse_args())
