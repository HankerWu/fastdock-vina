"""Compare common-reference scores for every frozen parent pose and self-reference for each first."""

import argparse
from pathlib import Path
import tempfile
from types import SimpleNamespace

from panel import read, write
from fixed import main as rescore, TERMS


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "results", "out", "parent"):
        p.add_argument(name, type=Path)
    p.add_argument("method")
    p.add_argument("candidate", type=Path)
    a = p.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        for item in ("runs.csv", "raw"):
            (d / item).symlink_to((a.results / item).resolve())
        rescore(
            SimpleNamespace(
                work=a.inputs,
                study="macro",
                output=d,
                parent=a.parent,
                candidate=a.candidate,
            )
        )
        rows = read(d / "fixed.csv")
        for parent, trial in zip(rows[::2], rows[1::2]):
            fields = ["key", "rank", "mode", "pose_sha256", *TERMS]
            if any(parent[k] != trial[k] for k in fields):
                raise ValueError(parent["key"] + ": frozen-coordinate scores diverged")
            trial["arm"] = a.method
        write(a.out, rows)


if __name__ == "__main__":
    main()
