"""Evaluate crystal-informed local starts outside the blind prediction panel."""

import argparse
from pathlib import Path

from panel import PANELS  # loads the matching graph/alignment helpers


def run(work, binary):
    from crystal_start import run as relax_starts

    pairs = {
        (query, next(iter(entries))) for entries in PANELS.values() for query in entries
    }
    relax_starts(work, binary, pairs, compare_mappings=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    args = p.parse_args()
    run(args.work.resolve(), args.binary.resolve())
