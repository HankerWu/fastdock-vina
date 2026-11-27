"""Run the independently registered high-water-mark structural experiment."""

import argparse
import csv
import re
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import audit_transfer
import prepare_inputs
import run_panel

SYSTEMS = {
    "1PMV": ("JNK3", "537", "A"),
    "2YXJ": ("BCL-XL", "N3C", "A"),
    "1T7T": ("AR", "DHT", "A"),
}
INSTANCES = {"1PMV": 501, "2YXJ": 1001, "1T7T": 300}
ARMS = ("upstream", "parent", "candidate", "hold")
SEEDS = (438013, 859337)

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "dock", "audit"))
    p.add_argument("work", type=Path)
    p.add_argument("--cpu", type=int, default=4)
    p.add_argument("--reverse", action="store_true")
    p.add_argument("--rescore", action="store_true")
    for arm in (*ARMS, "counter", "hold-counter"):
        p.add_argument("--" + arm, type=Path)
    args = p.parse_args()
    if args.action == "prepare":
        prepare_inputs.main(args.work, SYSTEMS)
    elif args.action == "dock":
        run_panel.run(args, SYSTEMS, SEEDS, ARMS)
    else:
        inputs = json.loads((args.work / "sources.json").read_text())
        eligible = {
            row["pdb"]: INSTANCES[row["pdb"]]
            for row in inputs
            if row["status"] == "ready"
        }
        audit_transfer.main(args.work, args.work, eligible, SYSTEMS)
        with (args.work / "cells.csv").open() as handle:
            for row in csv.DictReader(handle):
                if row["arm"] == "hold":
                    log = (args.work / "runs" / (row["key"] + ".counts")).read_text()
                    skipped = [
                        int(n) for n in re.findall(r"Skipped refills: (\d+)", log)
                    ]
                    assert len(skipped) == 4 and sum(skipped) == int(
                        row["skipped_refills"]
                    )
                    assert sum(n > 0 for n in skipped) == int(
                        row["trajectories_skipping_refills"]
                    )
