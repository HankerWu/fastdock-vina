"""Build preregistered MAPK14 structures using the frozen dry preparation."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "cumulative/external"))
import prepare as historical

CASES = [("1BL7", "SB4", "MAPK14"), ("1DI9", "MSQ", "MAPK14")]
historical.CASES = CASES

if __name__ == "__main__":
    historical.main()
