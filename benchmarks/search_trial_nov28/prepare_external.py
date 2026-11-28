"""Run the unchanged historical chemistry pipeline on the newly registered deposits."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "cumulative/external"))
import prepare

prepare.CASES = [
    ("4BJX", "73B", "BRD4"),
    ("4WIV", "3P2", "BRD4"),
    ("5U2E", "837", "BRD4"),
]
if __name__ == "__main__":
    prepare.main()
