"""Alternate before/after in-memory parser probes on the two FXa preparations."""

import argparse
import csv
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "old", "new"):
        parser.add_argument(name, type=lambda text: Path(text).resolve())
    args = parser.parse_args()
    rows = []
    for repeat in range(10):
        for geometry in ("original", "random_coords"):
            ligand = args.run / (geometry + ".pdbqt")
            order = ("new", "old") if repeat % 2 else ("old", "new")
            for arm in order:
                result = subprocess.check_output(
                    [getattr(args, arm), "time", "XS", "300", ligand], text=True
                )
                checksum, cpu, wall = result.split()
                rows.append(
                    dict(
                        round=repeat,
                        geometry=geometry,
                        arm=arm,
                        repeats=300,
                        pair_checksum=checksum,
                        cpu_s=cpu,
                        wall_s=wall,
                    )
                )
    with (args.run / "time.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
