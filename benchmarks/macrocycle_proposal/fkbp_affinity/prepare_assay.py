"""Prepare the Table 2 human FKBP12 Kd controls from deposited CCD graphs."""

import argparse
import os
import sys
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("repo", type=Path)
p.add_argument("work", type=Path)
p.add_argument("ids", nargs="*", default=["1FKJ", "6VCU"])
a = p.parse_args()
sys.path.insert(0, str(a.repo.resolve() / "benchmarks/pair_gradient/pose"))
import prepare

os.environ.setdefault("OPENMM_CPU_THREADS", "1")
for pdb in a.ids:
    prepare.CASES[pdb] = (
        "human FKBP12",
        {"1FKJ": "FK5", "6VCU": "R27", "1FKF": "FK5"}[pdb],
        "A",
    )
    prepare.prepare(a.work, pdb)
