"""Retain the fixed search prior for ranking previously native-refined scaffold pools."""

import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np
from ammonium_states import LIGANDS, RECEPTORS, SEEDS, xyzs
from coordination_probe import read_models, write_csv


def rank(runs, metadata):
    by_key = {}
    for row in csv.DictReader((runs / "candidates.csv").open()):
        by_key[
            (
                row["receptor"],
                row["ligand"],
                int(row["seed"]),
                row["slot"],
                int(row["model"]),
            )
        ] = row
    for path in sorted(runs.glob("*.pdbqt")):
        rec, lig, seed, slot = path.stem.split("_", maxsplit=3)
        sites = np.loadtxt(runs / "inputs" / f"{rec}_{lig}.anchor", ndmin=2)
        if sites.shape[1] != 4:
            raise ValueError("Only points are supported in this scaffold diagnostic")
        for m, lines in read_models(path.read_text().splitlines()).items():
            row = by_key[rec, lig, int(seed), slot, m]
            atoms = xyzs("\n".join(lines))
            d2 = np.sum(
                (np.array([atoms[int(s)] for s in sites[:, 0]]) - sites[:, 1:]) ** 2,
                axis=1,
            )
            # 4*(sqrt(1+d²)-1) is algebraically identical to the search's stable rational expression.
            row["prior"] = 4.0 * float(np.mean(np.sqrt(1.0 + d2) - 1.0))
    cells = []
    modes = []
    for rec, lig, seed in itertools.product(RECEPTORS, LIGANDS, SEEDS):
        for arm, second in (
            ("original_only", "second_original"),
            ("both_states", "second_reflected"),
        ):
            pool = [
                r
                for (q, l, s, slot, _), r in by_key.items()
                if (q, l, s) == (rec, lig, seed) and slot in ("common", second)
            ]

            def key(r, method):
                total = float(r["total"])
                score = {
                    "release": total,
                    "retain": total + r["prior"],
                    "geometry_only": r["prior"],
                }[method]
                return (score, r["slot"] != "common", int(r["model"]))

            raw_best = min(pool, key=lambda r: key(r, "release"))
            unbound = float(raw_best["intra"])
            for r in pool:
                modes.append(
                    {
                        "receptor": rec,
                        "ligand": lig,
                        "seed": seed,
                        "arm": arm,
                        "slot": r["slot"],
                        "model": int(r["model"]),
                        "native_raw": float(r["total"]),
                        "prior": r["prior"],
                        "native_plus_prior": float(r["total"]) + r["prior"],
                    }
                )
            for method in ("release", "retain", "geometry_only"):
                chosen = min(pool, key=lambda r: key(r, method))
                cells.append(
                    {
                        "receptor": rec,
                        "ligand": lig,
                        "seed": seed,
                        "arm": arm,
                        "method": method,
                        "native_affinity": (float(chosen["total"]) - unbound)
                        * float(chosen["flex_factor"]),
                        "rmsd": float(chosen["rmsd"]),
                        "sign": int(chosen["same_sign_as_bound"]),
                        "oracle": min(float(r["rmsd"]) for r in pool),
                        "oracle_signed": min(
                            (
                                float(r["rmsd"])
                                for r in pool
                                if int(r["same_sign_as_bound"])
                            ),
                            default=None,
                        ),
                        "slot": chosen["slot"],
                        "model": int(chosen["model"]),
                        "prior": chosen["prior"],
                        "common_unbound": unbound,
                    }
                )
    pairs = []
    for rec, seed, arm, method in itertools.product(
        RECEPTORS,
        SEEDS,
        ("original_only", "both_states"),
        ("release", "retain", "geometry_only"),
    ):
        group = {
            r["ligand"]: r["native_affinity"]
            for r in cells
            if (r["receptor"], r["seed"], r["arm"], r["method"])
            == (rec, seed, arm, method)
        }
        for a, b in itertools.combinations(
            sorted(LIGANDS, key=lambda l: metadata[l]["Ki_nM"]), 2
        ):
            if metadata[b]["Ki_nM"] / metadata[a]["Ki_nM"] < 3:
                continue
            pairs.append(
                {
                    "receptor": rec,
                    "seed": seed,
                    "arm": arm,
                    "method": method,
                    "stronger": a,
                    "weaker": b,
                    "correct": 0.5
                    if group[a] == group[b]
                    else int(group[a] < group[b]),
                }
            )
    for name, table in (
        ("rank_energies", modes),
        ("rank_cells", cells),
        ("rank_pairs", pairs),
    ):
        write_csv(runs / f"{name}.csv", table)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    args = parser.parse_args()
    rank(args.runs, json.loads((args.prepared / "metadata.json").read_text()))
