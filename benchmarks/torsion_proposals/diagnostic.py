"""Localize search versus fixed-geometry errors by relaxing crystal-derived starting poses."""

import argparse
import sys
from pathlib import Path

from paired import aligned_receivers


def main(repo, work, binary):
    sys.path.insert(0, str(repo.resolve() / "benchmarks/paired_crossdocking"))
    import crystal_start
    import run

    crystal_start.alignment = lambda a, b: run.alignment(*aligned_receivers(a, b))
    crystal_start.run(work, binary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for field in ("repo", "work", "binary"):
        p.add_argument(field, type=Path)
    a = p.parse_args()
    main(a.repo, a.work, a.binary.resolve())
