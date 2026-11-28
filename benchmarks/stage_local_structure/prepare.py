"""Prepare the prospectively selected ESR1 deposits under the historical inclusion rules."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "cumulative/external"))
import prepare as historical

CASES = [("1QKT", "EST", "ESR1"), ("1X7R", "GEN", "ESR1")]
historical.CASES = CASES

if __name__ == "__main__":
    historical.main()
