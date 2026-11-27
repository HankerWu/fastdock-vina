"""Diagnose local stability with crystallographic and torsion-projected starts."""

import argparse
import csv
import json
import sys
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "work", "parent", "fixed"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    w = a.work.resolve()
    sys.path.insert(0, str(a.repo.resolve() / "benchmarks/paired_crossdocking"))
    import crystal_start

    with (w / "cells.csv").open() as f:
        pairs = {
            (r["query"], r["receiver"])
            for r in csv.DictReader(f)
            if r["query"] in {"4W9X", "6VN8"}
        }
    codes = set().union(*pairs)
    meta = {q: json.loads((w / q / "preparation.json").read_text()) for q in codes}
    for engine in ("parent", "fixed"):
        d = w / "guided" / engine
        d.mkdir(parents=True, exist_ok=True)
        (d / "metadata.json").write_text(json.dumps(meta))
        for code in codes:
            if not (d / code).exists():
                (d / code).symlink_to(w / code, target_is_directory=True)
        crystal_start.run(d, getattr(a, engine).resolve(), pairs, compare_mappings=True)


if __name__ == "__main__":
    main()
