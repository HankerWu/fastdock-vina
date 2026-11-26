"""Convergence of the hydroxyl-rotation integral on the identity-unambiguous 4.79-fold edge."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from itertools import product, chain
from evaluate import thermal
from data import SEEDS, write
from quantum import binding_site, energy, ligand, merge, read_models, rotations


def measure(args):
    p, l, angles, rec, name, seed = args
    vals = []
    for a, pg, ep in angles:
        ec = energy(merge(pg, l), "water")
        vals.append(
            dict(
                receptor=rec, ligand=name, seed=seed, angle_deg=a, pocket=ep, complex=ec
            )
        )
    rows = []
    for stride in (8, 4, 2, 1):
        sel = vals[::stride]
        e = (
            thermal([r["complex"] for r in sel])
            - thermal([r["pocket"] for r in sel])
            - p
        )
        rows.append(
            dict(
                receptor=rec, ligand=name, seed=seed, step_deg=7.5 * stride, transfer=e
            )
        )
    return vals, rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("prepared", type=Path)
    p.add_argument("results", type=Path)
    p.add_argument("--jobs", type=int, default=4)
    a = p.parse_args()
    rec_angles = {}
    phis = [i * 7.5 for i in range(48)]
    for rec in ("1BTY", "3PTB"):
        ps = rotations(
            binding_site(a.prepared / (rec + ".pdbqt")),
            ("triad", "delta", "water"),
            phis,
        )
        rec_angles[rec] = [
            (phi, g, energy(g, "water")) for phi, g in zip(phis, ps, strict=True)
        ]
    jobs = []
    for rec, name, seed in product(rec_angles, ("NH2", "OCH3"), SEEDS):
        l = ligand(
            read_models(
                (a.prepared / f"{rec}_{name}_{seed}_shell.pdbqt")
                .read_text()
                .splitlines()
            )[1]
        )
        jobs.append((energy(l, "water"), l, rec_angles[rec], rec, name, seed))
    with ProcessPoolExecutor(a.jobs) as exe:
        results = list(exe.map(measure, jobs))
    write(
        a.results / "angular_raw.csv",
        list(chain.from_iterable(pair[0] for pair in results)),
    )
    write(
        a.results / "angular.csv",
        list(chain.from_iterable(pair[1] for pair in results)),
    )
