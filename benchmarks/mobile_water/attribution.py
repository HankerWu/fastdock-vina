"""Identify dry hits lost to the mobile term and the exact dry-winner/challenger crossover.

The crossover is an after-the-fact sensitivity diagnostic for the single site.
It is not a new prior for scoring other cases.
"""

import argparse
import math
from pathlib import Path

from mobile_water import PRESENT, RT
from replay import records, save


def blocked(row):
    return -math.expm1(-(float(row["mobile"]) - float(row["dry"])) / RT) / PRESENT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    directory = parser.parse_args().results
    groups = {}
    for r in records(directory / "poses.csv"):
        if r["target"] == "thrombin" and r["rmsd"]:
            groups.setdefault((r["receiver"], r["ligand"], r["seed"]), []).append(r)
    out = []
    for key, pool in groups.items():
        native = min(pool, key=lambda r: (float(r["dry"]), int(r["model"])))
        changed = min(pool, key=lambda r: (float(r["mobile"]), int(r["model"])))
        if float(native["rmsd"]) > 2 or float(changed["rmsd"]) <= 2:
            continue
        ratio = math.exp((float(native["dry"]) - float(changed["dry"])) / RT)
        crossover = (1 - ratio) / (blocked(native) - ratio * blocked(changed))
        assert 0 < crossover < PRESENT
        out.append(
            {
                "receiver": key[0],
                "ligand": key[1],
                "seed": key[2],
                "dry_winner": native["model"],
                "mobile_winner": changed["model"],
                "dry_gap": float(changed["dry"]) - float(native["dry"]),
                "penalty_dry_winner": float(native["mobile"]) - float(native["dry"]),
                "penalty_challenger": float(changed["mobile"]) - float(changed["dry"]),
                "pair_crossover_prior": crossover,
            }
        )
    save(directory / "lost_hits.csv", out)


if __name__ == "__main__":
    main()
