"""Locate the closest returned pair and recover its exact pre-refinement distances by slot."""

import csv
import itertools
import re
import sys
from pathlib import Path

import numpy as np

STAGE = re.compile(r"^\**STAGE (before|after) (\d+)\n(.*?)^END STAGE$", re.M | re.S)


def inspect(path):
    xyz = {}
    for stage, slot, block in STAGE.findall(path.read_text()):
        coordinates = block.split("PDBQT\n")[0].splitlines()[1:]
        xyz[stage, int(slot)] = np.array(
            [[float.fromhex(v) for v in row.split()] for row in coordinates]
        )

    def distance(stage, a, b):
        return np.sqrt(np.square(xyz[stage, a] - xyz[stage, b]).sum(axis=1).mean())

    pairs = list(
        itertools.combinations(sorted(k for stage, k in xyz if stage == "after"), 2)
    )
    if not pairs:
        raise ValueError(f"No eligible pose pair in {path}")
    a, b = min(pairs, key=lambda ij: distance("after", *ij))
    return [
        path.stem,
        len(pairs),
        a,
        b,
        distance("before", a, b),
        distance("after", a, b),
        min(distance("before", i, j) for i, j in pairs),
        sum(distance("after", i, j) < 1 for i, j in pairs),
    ]


if __name__ == "__main__":
    output = csv.writer(sys.stdout, lineterminator="\n")
    output.writerow(
        (
            "condition",
            "pairs",
            "slot_a",
            "slot_b",
            "same_pair_before_A",
            "minimum_after_A",
            "minimum_before_A",
            "pairs_below_1A_after",
        )
    )
    output.writerows(inspect(path) for path in sorted(Path(sys.argv[1]).glob("*.log")))
