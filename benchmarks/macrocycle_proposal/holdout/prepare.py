"""Screen the three structures registered for the persistent-ring confirmation."""

import argparse
from pathlib import Path
import sys

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("work", type=Path)
a = p.parse_args()
here = Path(__file__).resolve().parent
sys.path.insert(0, str(here.parent / "external"))
from prepare_panel import CASES, run

extra = {
    "5K0X": ("MERTK", "K0X", "A"),
    "3OWN": ("REN", "3OX", "A"),
    "5N1Z": ("BCL6", "8GQ", "A"),
}
CASES.update(extra)
for pdb in extra:
    run(here.parents[2], a.work.resolve(), pdb)
