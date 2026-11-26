"""Rank an existing native pose pool with a user-selected O--C--N torsion penalty."""

import argparse
import csv
import math
import sys
from pathlib import Path

from amide_rotamers import active_amide_sites, check_barrier, dihedral_degrees
from score_physical_pool import pdbqt_chemistry_record
from select_zinc_pose import parse_atoms, read_models, read_score_header


def rank(ligand, pool, barrier, geometry="distinguished"):
    check_barrier(barrier)
    tree = pdbqt_chemistry_record(ligand)
    sites = active_amide_sites(ligand, geometry)
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
        angles = [
            [dihedral_degrees((*site[:3], tip), coords) for tip in site[3:]]
            for site in sites
        ]
        penalty = barrier * sum(
            sum(math.sin(math.radians(a)) ** 2 for a in arms) / len(arms)
            for arms in angles
        )
        rows.append(
            {
                "model": number,
                "native": header["native_score"],
                "amide_penalty": penalty,
                "score": header["native_score"] + penalty,
                "dihedrals_deg": " ".join(
                    ",".join(f"{a:.4f}" for a in arms) for arms in angles
                ),
            }
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
    p.add_argument(
        "--geometry",
        choices=("distinguished", "all-arms"),
        default="distinguished",
        help="historical single substituent or mean sin² across all N-carbon arms of each movable non-imide site",
    )
    a = p.parse_args()
    rows = rank(a.ligand.read_text(), a.pool.read_text(), a.barrier, a.geometry)
    w = csv.DictWriter(sys.stdout, fieldnames=rows[0], lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
