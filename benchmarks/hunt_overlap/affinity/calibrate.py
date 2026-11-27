"""Budget interpolation using docking wall time alone."""

import argparse
from pathlib import Path

from assay import read


def cap(folder):
    def total(arm, budget):
        rows = read(folder / f"runs-912737-{budget}-{arm}.csv")
        assert len(rows) == 52
        return sum(float(r["seconds"]) for r in rows)

    reference = total("parent", 12000)
    low, high = total("candidate", 6000), total("candidate", 12000)
    if high <= low:
        # A flat calibration cannot resolve a higher budget at equal time.
        return 6000
    estimate = 6000 + 6000 * (reference - low) / (high - low)
    return int(min(12000, max(6000, 100 * round(estimate / 100))))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    print(cap(p.parse_args().folder))
