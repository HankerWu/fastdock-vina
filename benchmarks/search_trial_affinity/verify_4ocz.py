"""Assemble the three archived EPHX2 pools and the new experimental pool, then re-read the deposits."""

import argparse
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
EXT = HERE.parent / "cumulative/external"
sys.path.insert(0, str(HERE.parent / "geometry_comparison/upstream/four_versions"))
from run_panel import csvfile
from evaluate import read, summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for source, destination in [
            (EXT / "inputs.tar.xz", tmp / "fixtures"),
            (EXT / "results/raw.tar.xz", tmp / "old"),
            (HERE / "controls/control_raw.tar.xz", tmp / "new"),
        ]:
            with tarfile.open(source) as t:
                t.extractall(destination, filter="data")
        data = {}
        for arm, origin in (
            ("official", "official"),
            ("parent", "parent"),
            ("stable", "candidate"),
            ("candidate", "candidate"),
        ):
            folder = tmp / "runs" / arm
            folder.mkdir(parents=True)
            fresh = arm == "candidate"
            directory = HERE / "controls" if fresh else EXT / "results"
            raw = tmp / ("new/fresh/raw" if fresh else "old")
            (folder / "raw").mkdir()
            for filename in ("cells.csv", "ranks.csv"):
                selected = []
                for row in read(directory / filename):
                    original = row["key"]
                    if not original.endswith("-" + origin):
                        continue
                    row["key"] = original.rsplit("-", 1)[0] + "-" + arm
                    if filename == "cells.csv":
                        row["arm"] = arm
                        (folder / "raw" / (row["key"] + ".pdbqt")).symlink_to(
                            raw / (original + ".pdbqt")
                        )
                    selected.append(row)
                data.setdefault(filename, []).extend(selected)
                csvfile(folder / filename, selected)
        subprocess.run(
            [sys.executable, EXT / "audit.py", tmp / "fixtures", tmp / "runs"],
            stdout=(out / "independent.json").open("w"),
            check=True,
        )
        for filename, values in data.items():
            csvfile(out / filename, values)
        summaries(out)


if __name__ == "__main__":
    main()
