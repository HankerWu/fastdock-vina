"""Summarize top-1 physical checks and compare common-reference diagnostic energies."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from summary import seed_cohort


def dump(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def main(work):
    cells = {}
    for r in csv.DictReader((work / "cells.csv").open()):
        cells[tuple(r[z] for z in ("cap", "ligand", "receiver", "seed", "arm"))] = r
    q = defaultdict(list)
    keys = set()
    for p in sorted(work.glob("physical-validity*.csv")):
        for r in csv.DictReader(p.open()):
            k = tuple(r[z] for z in ("cap", "ligand", "receiver", "seed", "arm"))
            if k in keys:
                raise ValueError("duplicate validity row")
            keys.add(k)
            cohort = seed_cohort(r["seed"])
            q[(r["cap"], cohort, r["arm"])].append(
                (int(r["valid"]), float(cells[k]["top1"]))
            )
    if set(cells) != keys:
        raise ValueError("missing raw physical checks")
    rows = [
        dict(
            cap=cap,
            cohort=cohort,
            arm=arm,
            n=len(g),
            pb_valid=sum(v for v, _ in g),
            joint=sum(v and rms <= 2 for v, rms in g),
        )
        for (cap, cohort, arm), g in sorted(q.items())
    ]
    dump(work / "joint.csv", rows)
    rows = []
    for r in csv.DictReader((work / "crystal_starts.csv").open()):
        if r["arm"] != "projected":
            continue
        a, b = r["source"], r["receiver"]
        energies = [
            float(
                next(
                    line.split()[-1]
                    for line in (
                        work / "runs_12000" / f"{a}_{b}_{seed}_native.pdbqt"
                    ).open()
                    if line.startswith("REMARK INTER + INTRA:")
                )
            )
            for seed in (51913, 51941, 51973)
        ]
        rows.append(
            dict(
                ligand=a,
                receiver=b,
                local_rmsd=r["rmsd"],
                local_energy=r["inter_plus_intra"],
                best_global_energy=min(energies),
                local_minus_best=float(r["inter_plus_intra"]) - min(energies),
            )
        )
    dump(work / "projected_gaps.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    main(p.parse_args().work)
