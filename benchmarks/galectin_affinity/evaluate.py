"""Audit raw native scores and the common pyranose core at every returned rank."""

import argparse
import csv
import itertools
import math
from pathlib import Path
from statistics import mean, median

import numpy as np
from dock import ARMS, SEEDS
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from scipy.stats import spearmanr
from sugar import DONORS, crystal_core, ligand_core, xyz

TARGETS = {"Gal1": ("1GZW", "4Q26"), "Gal3": ("1KJL", "3T1L")}
RT = 0.00198720425864083 * 298.15


def write(path, rows):
    with path.open("w") as outf:
        writer = csv.DictWriter(outf, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def poses(work):
    cores = {p: crystal_core(work, p) for p in DONORS}
    serials = {i: ligand_core(work, i) for i in range(5, 18)}
    rows = []
    details = []
    for pdb, compound, seed, arm in itertools.product(
        DONORS, range(5, 18), SEEDS, ARMS
    ):
        name = f"{pdb}_{compound}_{seed}_{arm}"
        text = (work / "pools" / f"{name}.pdbqt").read_text()
        cell = []
        for mode, chunk in enumerate(text.split("ENDMDL")[:-1], 1):
            lines = chunk.strip().splitlines()
            score = float(
                next(r.split()[3] for r in lines if r.startswith("REMARK VINA RESULT:"))
            )
            atoms = {
                int(s[6:11]): xyz(s)
                for s in lines
                if s.startswith(("ATOM  ", "HETATM"))
            }
            delta = np.array([atoms[a] for a in serials[compound]]) - cores[pdb]
            rmsd = float(np.sqrt(np.mean(np.sum(delta**2, axis=1))))
            cell.append((score, rmsd))
            details.append(
                {
                    "pdb": pdb,
                    "compound": compound,
                    "seed": seed,
                    "arm": arm,
                    "mode": mode,
                    "score": score,
                    "core_rmsd": rmsd,
                }
            )
        supported = [(s, r) for s, r in cell if r <= 2.0]
        rows.append(
            {
                "pdb": pdb,
                "compound": compound,
                "seed": seed,
                "arm": arm,
                "modes": len(cell),
                "top": cell[0][0],
                "core_conditioned": min(s for s, r in supported)
                if supported
                else cell[0][0],
                "top_rmsd": cell[0][1],
                "best_rmsd": min(r for s, r in cell),
                "any_core": bool(supported),
            }
        )
    return rows, details


def aggregate(states, scheme):
    if scheme == "mean":
        return mean(states)
    if scheme == "exponential":
        m = min(states)
        return m - RT * math.log(mean(math.exp((m - s) / RT) for s in states))
    return states[int(scheme)]


def evaluate(work, out):
    out.mkdir(exist_ok=True, parents=True)
    scores, details = poses(work)
    write(out / "cells.csv", scores)
    write(out / "poses.csv", details)
    bycell = {(r["pdb"], r["compound"], r["arm"], r["seed"]): r for r in scores}
    panel = {int(r["compound"]): r for r in csv.DictReader((work / "panel.csv").open())}
    for row in panel.values():
        row["formula"] = rdMolDescriptors.CalcMolFormula(
            Chem.MolFromSmiles(row["smiles"])
        )
    summaries = []
    pairs = []
    geometries = []
    bootstrap = []
    for target, pdbs in TARGETS.items():
        kd = {i: float(row[f"{target}_Kd_uM"]) for i, row in panel.items()}
        for arm in ARMS:
            for pdb in pdbs:
                g = [r for r in scores if r["arm"] == arm and r["pdb"] == pdb]
                geometries.append(
                    {
                        "target": target,
                        "pdb": pdb,
                        "arm": arm,
                        "n": len(g),
                        "top_core": sum(r["top_rmsd"] <= 2 for r in g),
                        "pool_core": sum(r["any_core"] for r in g),
                    }
                )
        for cohort, seeds in [
            ("all", SEEDS),
            ("first", SEEDS[:2]),
            ("second", SEEDS[2:]),
        ]:
            for field, scheme in itertools.product(
                ("top", "core_conditioned"), ("mean", "exponential", "0", "1")
            ):
                energies = {
                    arm: {
                        i: aggregate(
                            [
                                median(bycell[pdb, i, arm, s][field] for s in seeds)
                                for pdb in pdbs
                            ],
                            scheme,
                        )
                        for i in kd
                    }
                    for arm in ARMS
                }
                energies["larger_heavy"] = {
                    i: -int(panel[i]["heavy_atoms"]) for i in kd
                }
                for subset in ("all", "triazoles", "same_formula"):
                    eligible = [
                        (i, j)
                        for i, j in itertools.permutations(kd, 2)
                        if kd[j] >= 3 * kd[i]
                        and (subset == "all" or 17 not in (i, j))
                        and (
                            subset != "same_formula"
                            or panel[i]["formula"] == panel[j]["formula"]
                        )
                    ]
                    results = {}
                    for arm, E in energies.items():
                        gaps = [E[j] - E[i] for i, j in eligible]
                        outcomes = np.asarray(
                            [1 if z > 1e-9 else 0 if z < -1e-9 else 0.5 for z in gaps]
                        )
                        results[arm] = outcomes
                        if subset == "same_formula":
                            rho = ""
                        else:
                            ids = [i for i in kd if subset == "all" or i != 17]
                            rho = float(
                                spearmanr(
                                    [kd[i] for i in ids], [E[i] for i in ids]
                                ).statistic
                            )
                        summaries.append(
                            {
                                "target": target,
                                "cohort": cohort,
                                "field": field,
                                "scheme": scheme,
                                "subset": subset,
                                "arm": arm,
                                "n": len(eligible),
                                "correct": float(outcomes.sum()),
                                "spearman": rho,
                            }
                        )
                        if (
                            cohort == "all"
                            and field == "top"
                            and scheme == "mean"
                            and subset == "all"
                        ):
                            for (i, j), gap, c in zip(eligible, gaps, outcomes):
                                pairs.append(
                                    {
                                        "target": target,
                                        "strong": i,
                                        "weak": j,
                                        "strong_kd": kd[i],
                                        "weak_kd": kd[j],
                                        "exp_gap": RT * math.log(kd[j] / kd[i]),
                                        "arm": arm,
                                        "score_strong": E[i],
                                        "score_weak": E[j],
                                        "score_gap": gap,
                                        "correct": float(c),
                                        "same_formula": panel[i]["formula"]
                                        == panel[j]["formula"],
                                    }
                                )
                    if (
                        field == "top"
                        and scheme == "mean"
                        and cohort == "all"
                        and subset in ("all", "triazoles")
                    ):
                        delta = results["balanced"] - results["native"]
                        ids = sorted(kd) if subset == "all" else sorted(set(kd) - {17})
                        rng = np.random.default_rng(20261126)
                        draws = []
                        for _ in range(4000):
                            w = np.bincount(rng.choice(ids, len(ids)), minlength=18)
                            weights = np.array([w[i] * w[j] for i, j in eligible])
                            if weights.sum():
                                draws.append(float(delta @ weights / weights.sum()))
                        lo, hi = np.quantile(draws, [0.025, 0.975])
                        bootstrap.append(
                            {
                                "target": target,
                                "subset": subset,
                                "delta": float(delta.mean()),
                                "lo": float(lo),
                                "hi": float(hi),
                                "draws": len(draws),
                            }
                        )
    write(out / "summary.csv", summaries)
    write(out / "pairs.csv", pairs)
    write(out / "core_counts.csv", geometries)
    write(out / "paired_bootstrap.csv", bootstrap)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    evaluate(args.work, args.out)
