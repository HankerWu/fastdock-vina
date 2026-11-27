"""Score a native output PDBQT against doubly anchored receiver waters."""

import argparse
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from mobile_water import Site
from replay import cross_sites, models, pdb_lines, save, xyz


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receptor-pdb", type=Path, required=True)
    parser.add_argument("--chain", required=True, help="receiving protein chain ID")
    parser.add_argument(
        "--center", type=float, nargs=3, required=True, metavar=("X", "Y", "Z")
    )
    parser.add_argument(
        "--pool",
        type=Path,
        required=True,
        help="single-ligand Vina output with MODEL records",
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    rows = pdb_lines(args.receptor_pdb)
    coords = [xyz(l) for l in rows if l.startswith("ATOM") and l[21] == args.chain]
    if not coords:
        raise ValueError(f"No protein atoms for chain {args.chain!r}")
    tree = cKDTree(coords)
    sites = cross_sites(rows, tree, np.array(args.center), args.chain)
    prep_seconds, quadratures = {}, {}
    for order in (11, 17):
        then = time.perf_counter()
        quadratures[order] = [
            Site.prepare(center, b, tree, order) for _, center, b in sites
        ]
        prep_seconds[order] = time.perf_counter() - then
    scored = []
    for model, dry, ligand in models(args.pool):
        energies, costs = {}, {}
        for order, grid in quadratures.items():
            then = time.perf_counter()
            energies[order] = [s.evaluate(ligand) for s in grid]
            costs[order] = time.perf_counter() - then
        fixed = sum(v[0] for v in energies[11])
        scored.append(
            {
                "model": model,
                "dry": dry,
                "point": dry + fixed,
                "mobile": dry + sum(v[1] for v in energies[11]),
                "order17": dry + sum(v[1] for v in energies[17]),
                "n_sites": len(sites),
                "seconds_order11": costs[11],
                "seconds_order17": costs[17],
            }
        )
    if not scored:
        raise ValueError("No MODEL records")
    save(args.out, scored)
    print(
        f"sites={len(sites)} receiver_order11_s={prep_seconds[11]:.6f} "
        f"receiver_order17_s={prep_seconds[17]:.6f} total_wall_s={time.perf_counter() - started:.6f}"
    )


if __name__ == "__main__":
    main()
