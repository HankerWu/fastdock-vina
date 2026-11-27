"""Independently reconstruct deposited chemical graphs and validate local reference poses."""

import argparse
import csv
import json
from pathlib import Path

from backtracking_panel import CASES


def check(work, parent):
    import audit
    from crystal_start import run as relax

    checks = audit.main(work, CASES)
    (work / "audit.json").write_text(json.dumps(checks, indent=2) + "\n")
    if parent:
        with (work / "reference_geometry.csv").open() as f:
            pairs = {(row["query"], row["receiver"]) for row in csv.DictReader(f)}
        meta = {
            pdb: json.loads((work / pdb / "preparation.json").read_text())
            for pdb in CASES
        }
        (work / "metadata.json").write_text(json.dumps(meta))
        relax(work, parent.resolve(), pairs, compare_mappings=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--parent", type=Path)
    args = p.parse_args()
    check(args.work, args.parent)
