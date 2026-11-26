"""Shared inputs and table format for the matched amidinium panels."""

import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
RT = 0.00198720425864 * (273.15 + 25.1)
RT_FXA = 0.00198720425864 * 298.15
SEEDS = (271828, 314159, 161803, 141421)
LIGANDS = ("H", "NH2", "CH3", "OCH3", "CONH2")


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
