"""Separate each solute's water and binding-site legs before pairing substitutions."""

import argparse
import itertools
from collections import defaultdict
from pathlib import Path
from statistics import median

import numpy as np
from scipy.special import logsumexp

from data import RT, read, write

HERE = Path(__file__).resolve().parent
KEY = ("receptor", "seed", "arm", "ligand")


def thermal(values):
    e = np.asarray(values, dtype=float)
    return float(-RT * (logsumexp(-e / RT) - np.log(len(e))))


def summarize(raw, obc):
    groups = defaultdict(list)
    for row in raw:
        groups[
            tuple(row[k] for k in (*KEY, "fragment", "histidine", "solvent"))
        ].append(row)
    cells = {}
    for allkeys, rs in groups.items():
        key = allkeys[:4]
        cell = cells.setdefault(key, dict(zip(KEY, key, strict=True)))
        group, hist, solv = allkeys[4:]
        rs.sort(key=lambda r: int(r["angle_deg"]))
        assert [int(r["angle_deg"]) for r in rs] == (
            [0] if group == "acetate" else list(range(0, 360, 60))
        )
        r = rs[0]
        cell["vina"] = float(r["native"])
        cell["water"] = float(r["ligand_water"]) - float(r["ligand_vacuum"])
        label = "_".join(allkeys[4:])
        pocket = [float(x["pocket"]) for x in rs]
        complexes = [float(x["complex"]) for x in rs]
        transfer = thermal(complexes) - thermal(pocket) - float(r["ligand_water"])
        cell[label] = transfer
        if group == "triad" and hist == "delta" and solv == "water":
            cell["primary_bound"] = transfer + cell["water"]
            cell["primary_fixed0"] = complexes[0] - pocket[0] - float(r["ligand_water"])
            cell["primary_best_bound_angle"] = int(
                rs[np.argmin(complexes)]["angle_deg"]
            )
            cell["primary_best_apo_angle"] = int(rs[np.argmin(pocket)]["angle_deg"])
    for row in obc:
        key = tuple(row[k] for k in KEY)
        cells[key]["gbsa_full_protein"] = float(row["polar_transfer"])
    return list(cells.values())


def differences(rows):
    series = read(HERE.parent / "directional_ion_pair/results/series.csv")
    measurements = {x["ligand"]: float(x["kd_uM"]) for x in series}
    byblock = defaultdict(dict)
    for r in rows:
        byblock[r["arm"], r["receptor"], r["seed"]][r["ligand"]] = r
    rows_out = []
    for (arm, receptor, seed), block in byblock.items():
        for a, b in itertools.combinations(sorted(block, key=measurements.get), 2):
            ratio = measurements[b] / measurements[a]
            result = dict(
                arm=arm,
                receptor=receptor,
                seed=seed,
                strong=a,
                weak=b,
                ratio=ratio,
                observed=RT * np.log(ratio),
            )
            for key in block[a]:
                if key not in KEY and not key.endswith("angle"):
                    result[key] = block[b][key] - block[a][key]
            rows_out.append(result)
    stats = []
    for arm in ("native", "shell"):
        for rec in ("1BTY", "3PTB", "both"):
            block = [
                r
                for r in rows_out
                if r["arm"] == arm and (rec == "both" or r["receptor"] == rec)
            ]
            for strong, weak in sorted({(r["strong"], r["weak"]) for r in block}):
                group = [r for r in block if (r["strong"], r["weak"]) == (strong, weak)]
                s = dict(
                    arm=arm,
                    receptor=rec,
                    strong=strong,
                    weak=weak,
                    matched_cells=len(group),
                    ratio=group[0]["ratio"],
                    observed=group[0]["observed"],
                )
                for key in group[0]:
                    if key not in (
                        "arm",
                        "receptor",
                        "seed",
                        "strong",
                        "weak",
                        "ratio",
                        "observed",
                    ):
                        s[key] = median([r[key] for r in group])
                stats.append(s)
    return rows_out, stats


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results", type=Path)
    a = p.parse_args()
    cells = summarize(
        read(a.results / "energies.csv"), read(a.results / "classical.csv")
    )
    write(a.results / "cells.csv", cells)
    pairs, summary = differences(cells)
    write(a.results / "paired.csv", pairs)
    write(a.results / "summary.csv", summary)
