"""Rank an existing native pose pool with a user-selected O--C--N torsion penalty."""

import argparse
import csv
import math
from pathlib import Path
import sys

from amide_rotamers import amides, dihedral_degrees
from score_physical_pool import pdbqt_chemistry_record
from select_zinc_pose import parse_atoms, read_models, read_score_header


def rank(ligand, pool, barrier):
    if not math.isfinite(barrier) or barrier < 0:
        raise ValueError("Barrier must be a finite, nonnegative energy in kcal/mol")
    tree = pdbqt_chemistry_record(ligand)
    active = {
        frozenset(map(int, l.split()[1:]))
        for l in ligand.splitlines()
        if l.startswith("BRANCH ")
    }
    sites = [s for s in amides(ligand) if frozenset(s[1:3]) in active]
    rows = []
    for number, lines in read_models(pool.splitlines()).items():
        if pdbqt_chemistry_record("\n".join(lines)) != tree:
            raise ValueError(
                f"MODEL {number} differs from the input graph, torsion tree or atom records"
            )
        header = read_score_header(number, lines)
        if header["pro_score"] is not None:
            raise ValueError("Expected a native Vina pose pool")
        coords = {i: a.xyz for i, a in parse_atoms(lines).items()}
        angles = [dihedral_degrees(s, coords) for s in sites]
        penalty = barrier * sum(math.sin(math.radians(a)) ** 2 for a in angles)
        rows.append(
            dict(
                model=number,
                native=header["native_score"],
                amide_penalty=penalty,
                score=header["native_score"] + penalty,
                dihedrals_deg=" ".join(f"{a:.4f}" for a in angles),
            )
        )
    if not rows:
        raise ValueError("Empty pose pool")
    return sorted(rows, key=lambda r: (r["score"], r["model"]))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("ligand", type=Path)
    p.add_argument("pool", type=Path)
    p.add_argument(
        "--barrier",
        type=float,
        required=True,
        help="energy at a 90-degree twist in kcal/mol",
    )
    a = p.parse_args()
    rows = rank(a.ligand.read_text(), a.pool.read_text(), a.barrier)
    w = csv.DictWriter(sys.stdout, fieldnames=rows[0], lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
