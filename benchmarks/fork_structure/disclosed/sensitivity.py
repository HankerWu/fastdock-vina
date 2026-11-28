"""Reconstruct the post hoc 32-cell reciprocal comparison, eight equally weighted target strata."""

import argparse
from collections import Counter
from pathlib import Path
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import evidence

SOURCES = ("fixed", "external", "disclosed/ESR1", "disclosed/BRD4")
FIELDS = ("target", "query", "receiver", "seed")


def collect(out):
    origins, cells, keys = [], [], set()
    for dataset in SOURCES:
        folder = HERE.parent / dataset
        for row in evidence.evaluate.read(folder / "cells.csv"):
            if row["cohort"] != "reciprocal":
                continue
            key = tuple(row[name] for name in ("arm", *FIELDS))
            assert key not in keys and row["key"] not in {r["key"] for r in cells}
            keys.add(key)
            origins.append(dict(key=row["key"], source=dataset + "/cells.csv"))
            cells.append(row)
    assert len(cells) == 128
    ref = None
    for arm in evidence.engine.ARMS:
        rows = [r for r in cells if r["arm"] == arm]
        assert len(rows) == 32
        counts = Counter(r["target"] for r in rows)
        assert len(counts) == 8 and set(counts.values()) == {4}
        assert set(r["seed"] for r in rows) == {"85231", "85237"}
        identities = {tuple(r[k] for k in FIELDS) for r in rows}
        assert ref is None or identities == ref
        ref = identities
    out.mkdir(parents=True, exist_ok=True)

    def ordered(row):
        return row["key"]

    evidence.engine.csvfile(out / "cells.csv", sorted(cells, key=ordered))
    evidence.engine.csvfile(out / "source_index.csv", sorted(origins, key=ordered))
    evidence.evaluate.summaries(out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    output = HERE / "pooled_reciprocal"
    if args.write:
        collect(output)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp)
            collect(dest)
            for name in (
                "cells.csv",
                "source_index.csv",
                "summary.csv",
                "paired.csv",
                "paired_summary.csv",
            ):
                assert (dest / name).read_bytes() == (output / name).read_bytes(), name
        print(
            "128 outputs, 32 unique recipient/query/seed cells per arm and eight equally weighted target strata verified."
        )
