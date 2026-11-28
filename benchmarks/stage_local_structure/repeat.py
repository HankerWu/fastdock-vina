"""Verify independent ESR1 reruns by coordinate hashes, stage counts and rank table bytes."""

import argparse
import json
from pathlib import Path
from external import collect, evaluate


def main(fixtures, work, primary, out):
    collect(fixtures, work, out)
    columns = (
        "status",
        "pool_sha256",
        "modes",
        "hunt_evals",
        "second_evals",
        "explicit_evals",
        "explicit_calls",
        "late_prescore",
    )

    def by_key(folder):
        return {
            r["key"]: tuple(r[k] for k in columns)
            for r in evaluate.read(folder / "cells.csv")
        }

    original = by_key(primary)
    if (
        by_key(out) != original
        or (out / "ranks.csv").read_bytes() != (primary / "ranks.csv").read_bytes()
    ):
        raise ValueError("repeated search or ranked geometry diverges")
    print(
        json.dumps(
            dict(
                identical_cells=len(original),
                identical_ranks=len(evaluate.read(out / "ranks.csv")),
                independent_geometry=json.loads(
                    (out / "verification.json").read_text()
                ),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("fixtures", "work", "primary", "out"):
        parser.add_argument(name, type=Path)
    main(**vars(parser.parse_args()))
