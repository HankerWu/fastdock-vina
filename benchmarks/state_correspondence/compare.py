"""Report paired two-state results relative to an exact declared-graph control."""

import argparse
import csv
import hashlib
import itertools
from pathlib import Path
from statistics import median


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def save(path, table):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(table[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(table)


def compare(control, fallback, released, out):
    out.mkdir(exist_ok=True, parents=True)
    summary = []
    output_hashes = []
    for scheme, folder in (
        ("declared", control),
        ("fallback", fallback),
        ("released_neighborhood", released),
    ):
        cells = rows(folder / "rank_cells.csv")
        pairs = rows(folder / "rank_pairs.csv")
        for arm, method in itertools.product(
            ("original_only", "both_states"), ("release", "retain", "geometry_only")
        ):
            c = [r for r in cells if (r["arm"], r["method"]) == (arm, method)]
            p = [r for r in pairs if (r["arm"], r["method"]) == (arm, method)]
            if len(c) != 20 or len(p) != 28:
                raise ValueError("Incorrect denominator")
            summary.append(
                {
                    "search": scheme,
                    "inputs": arm,
                    "selection": method,
                    "n": len(c),
                    "top2A": sum(float(r["rmsd"]) <= 2 for r in c),
                    "joint_top": sum(
                        float(r["rmsd"]) <= 2 and int(r["sign"]) for r in c
                    ),
                    "oracle": sum(float(r["oracle"]) <= 2 for r in c),
                    "joint_oracle": sum(
                        bool(r["oracle_signed"]) and float(r["oracle_signed"]) <= 2
                        for r in c
                    ),
                    "median_top_rmsd": median(float(r["rmsd"]) for r in c),
                    "Ki_pairs": len(p),
                    "Ki_credit": sum(float(r["correct"]) for r in p),
                }
            )
        for file in sorted(folder.glob("*.pdbqt")):
            digest = hashlib.sha256(file.read_bytes()).hexdigest()
            same = file.read_bytes() == (control / file.name).read_bytes()
            output_hashes.append(
                {
                    "search": scheme,
                    "filename": file.name,
                    "sha256": digest,
                    "identical_to_control": int(same),
                }
            )
            if scheme == "fallback" and not same:
                raise ValueError("Unexpected mismatch in fallback control")
    save(out / "summary.csv", summary)
    save(out / "pool_hashes.csv", output_hashes)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("control", "fallback", "released", "out"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    compare(a.control, a.fallback, a.released, a.out)
