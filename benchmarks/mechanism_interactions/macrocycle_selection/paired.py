"""Emit fixed-minus-parent differences in the order of the parent run manifest."""
import argparse
import csv
from pathlib import Path


FIELDS = (
    "modes", "severe_all", "inversions", "saved_top1_A", "saved_best5_A", "oracle_A",
    "saved_regret_A", "hunt_evals", "save_evals", "explicit_evals", "cpu_s", "wall_s",
)
CONDITIONS = ("pdb", "scoring", "stage", "seed", "trajectories", "cap_per_trajectory", "requested_modes")


def read_cells(path):
    with path.open(newline="") as stream:
        records = list(csv.DictReader(stream))
    by_key = {row["cell"]: row for row in records}
    if not records or len(by_key) != len(records):
        raise ValueError(f"Missing cells or duplicate cell identifiers in {path}")
    return by_key


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("parent_cells", "fixed_cells", "output"):
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    before, after = read_cells(args.parent_cells), read_cells(args.fixed_cells)
    if before.keys() != after.keys():
        raise ValueError("Parent and fixed tables have different cell identifiers")
    rows = []
    for cell, a in before.items():
        b = after[cell]
        if any(a[k] != b[k] for k in CONDITIONS):
            raise ValueError(f"{cell}: conditions changed")
        # Keep unavailable structural metrics missing, including when both runs return no poses.
        differences = [round(float(b[k]) - float(a[k]), 6) if a[k] and b[k] else "" for k in FIELDS]
        rows.append([cell, *differences])
    with args.output.open("w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["cell", *("delta_" + k for k in FIELDS)])
        writer.writerows(rows)


if __name__ == "__main__":
    main()
