"""Measure native-score ranking and symmetry-corrected receptor-frame ligand RMSD."""

import argparse
import csv
import itertools
import json
import math
from pathlib import Path
from statistics import median

import numpy as np
from run import ARMS, SEEDS, STATES

KD_UM = {"4W52": 175, "4W53": 102, "4W54": 68, "4W55": 18}  # Merski et al., Table 1


def write(path, rows):
    with path.open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def read_poses(path):
    for block in path.read_text().split("MODEL ")[1:]:
        if not block.strip().endswith("ENDMDL"):
            raise ValueError(f"Incomplete pose in {path}")
        lines = block.splitlines()
        rank = int(lines[0])
        (score,) = [
            float(x.split()[3]) for x in lines if x.startswith("REMARK VINA RESULT:")
        ]
        positions = {
            int(x[6:11]): [float(x[k : k + 8]) for k in (30, 38, 46)]
            for x in lines
            if x.startswith(("ATOM  ", "HETATM"))
        }
        yield rank, score, positions


def evaluate(work, out):
    out.mkdir(exist_ok=True, parents=True)
    prep = work / "prepared"
    manifest = json.loads((prep / "manifest.json").read_text())
    cells = []
    modes = []
    pairs = []
    for state, seed, arm in itertools.product(STATES, SEEDS, ARMS):
        group = []
        for lig, m in manifest.items():
            base = {"receptor_state": state, "seed": seed, "arm": arm, "ligand": lig}
            pool = []
            for rank, score, position in read_poses(
                work / "runs" / f"{state}_{lig}_{seed}_{arm}.pdbqt"
            ):
                xyz = np.array([position[s] for s in m["serials"]])
                by_alt = {
                    alt: float(
                        np.sqrt(
                            np.square(np.array(data["mapped"]) - xyz)
                            .sum(axis=-1)
                            .mean(axis=-1)
                        ).min()
                    )
                    for alt, data in m["refs"].items()
                }
                major = max(m["refs"], key=lambda a: m["refs"][a]["occupancy"])
                row = dict(
                    base,
                    rank=rank,
                    score=score,
                    rmsd=min(by_alt.values()),
                    major_rmsd=by_alt[major],
                    closest_alt=min(by_alt, key=by_alt.get),
                    all_alt_rmsd=json.dumps(by_alt, sort_keys=True),
                )
                modes.append(row)
                pool.append(row)
            if not 1 <= len(pool) <= 9 or [v["rank"] for v in pool] != list(
                range(1, len(pool) + 1)
            ):
                raise ValueError(f"Incomplete set at {base}")
            group.append(
                dict(
                    base,
                    top_score=pool[0]["score"],
                    top_rmsd=pool[0]["rmsd"],
                    top_major_rmsd=pool[0]["major_rmsd"],
                    top_alt=pool[0]["closest_alt"],
                    best_rmsd=min(v["rmsd"] for v in pool),
                    n_modes=len(pool),
                )
            )
        cells += group
        by_lig = {r["ligand"]: r for r in group}
        for a, b in itertools.combinations(KD_UM, 2):
            strong, weak = sorted((a, b), key=KD_UM.get)
            ratio = KD_UM[weak] / KD_UM[strong]
            ratio_from_dg = math.exp(
                (manifest[weak]["dg"] - manifest[strong]["dg"])
                / (0.00198720425864 * 302)
            )
            if (ratio >= 3) != (ratio_from_dg >= 3):
                raise ValueError("Printed Kd and free energies disagree on inclusion")
            if ratio < 3:
                continue
            delta = by_lig[weak]["top_score"] - by_lig[strong]["top_score"]
            pairs.append(
                {
                    "receptor_state": state,
                    "seed": seed,
                    "arm": arm,
                    "stronger": strong,
                    "weaker": weak,
                    "kd_ratio": ratio,
                    "score_margin": delta,
                    "credit": 1 if delta > 0 else 0.5 if delta == 0 else 0,
                }
            )
    summary = []
    for state, arm in itertools.product((*STATES, "all"), ARMS):

        def select(x, arm=arm, state=state):
            return x["arm"] == arm and (state == "all" or state == x["receptor_state"])

        c = list(filter(select, cells))
        p = list(filter(select, pairs))
        summary.append(
            {
                "receptor_state": state,
                "arm": arm,
                "cells": len(c),
                "top2A": sum(x["top_rmsd"] <= 2 for x in c),
                "major_top2A": sum(x["top_major_rmsd"] <= 2 for x in c),
                "pool2A": sum(x["best_rmsd"] <= 2 for x in c),
                "top_median_rmsd": median(x["top_rmsd"] for x in c),
                "pairs": len(p),
                "pair_credit": sum(x["credit"] for x in p),
            }
        )
    for name, rows in (
        ("poses", modes),
        ("cells", cells),
        ("pairs", pairs),
        ("summary", summary),
    ):
        write(out / f"{name}.csv", rows)
    print(f"{len(modes)} poses, {len(cells)} cells")
    print((out / "summary.csv").read_text())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    evaluate(a.work, a.out)
