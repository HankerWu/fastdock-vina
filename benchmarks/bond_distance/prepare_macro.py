"""Prepare two FKBP12 macrocycles with the archived, crystal-independent ligand recipe."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pair_gradient/pose"))
import prepare

if __name__ == "__main__":
    for pdb, component in (("1FKF", "FK5"), ("1FKB", "RAP")):
        prepare.CASES[pdb] = ("FKBP12", component, "A")
        prepare.prepare(Path(sys.argv[1]), pdb)
