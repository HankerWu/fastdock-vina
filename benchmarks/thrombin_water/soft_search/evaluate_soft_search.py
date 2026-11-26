"""Compare matched hard and soft searches on both seed panels; recompute from raw coordinates."""

import argparse
import csv
import itertools
import math
import sys
from pathlib import Path
from statistics import median

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from independent_audit import (
    get_serials,
    move_points,
    positions,
    read_csv,
    refs_from_ccd,
)

R_KJ = 0.00831446261815324
TEMP = 298.15
SEEDS = {"first": ("83519", "194027"), "confirmation": ("88091", "302123")}
RECEIVERS = ("3D49", "2ZFF", "2ZC9")


def dump(path, rows):
    with path.open("w") as h:
        w = csv.DictWriter(h, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)


def evaluate(phase, hard, soft, confirm, out, series):
    water = [[12.575, -7.936, 24.469]]  # independent 3D49 H:2006 deposited oxygen
    cells = []
    for r in RECEIVERS:
        site = move_points(water, "3D49", r)[0]
        for lig, exp in series.items():
            top = Path("prepared", lig + ".pdbqt").read_text().splitlines()
            serials = get_serials(top)
            p = exp["PDB"]
            refs = refs_from_ccd(p, exp["CCD"], top, r) if p else []
            for seed, arm in itertools.product(SEEDS[phase], ("hard", "soft")):
                if phase == "first":
                    path = (hard if arm == "hard" else soft) / f"{r}_{lig}_{seed}.pdbqt"
                else:
                    path = confirm / f"{r}_{lig}_{seed}_{arm}.pdbqt"
                pool = []
                for model in path.read_text().split("ENDMDL"):
                    if not model.strip():
                        continue
                    ls = model.splitlines()
                    mapping = {
                        int(l[6:11]): positions(l)
                        for l in ls
                        if l.startswith(("ATOM", "HETATM"))
                    }
                    xyz = np.array([mapping[i] for i in serials])
                    rms = min(
                        (
                            np.linalg.norm(xyz - ref) / math.sqrt(len(xyz))
                            for ref in refs
                        ),
                        default=None,
                    )
                    dist = min(math.dist(v, site) for v in xyz)
                    t = min(1, max(0, (dist - 2.0) / 0.6))
                    a = t * t * (3 - 2 * t)
                    native = float(
                        next(
                            l for l in ls if l.startswith("REMARK VINA RESULT")
                        ).split()[3]
                    )
                    wet = native - R_KJ * TEMP / 4.184 * math.log(0.25 + 0.75 * a)
                    pool.append((native, wet, rms))
                for k, score in enumerate(("native", "exclusion")):
                    winner = min(pool, key=lambda x: x[k])
                    cells.append(
                        {
                            "phase": phase,
                            "receiver": r,
                            "ligand": lig,
                            "seed": seed,
                            "arm": arm,
                            "score": score,
                            "energy": winner[k],
                            "rmsd": winner[2],
                            "best_rmsd": min(v[2] for v in pool) if p else None,
                            "kind": "cognate"
                            if p == r
                            else "cross"
                            if p
                            else "no_crystal",
                        }
                    )
    return cells


def summarize(cells, series, out):
    summary = []
    paired = []
    pairs = []
    omit = []
    for phase, arm, score in itertools.product(
        ("first", "confirmation", "all"), ("hard", "soft"), ("native", "exclusion")
    ):
        sel = [
            r
            for r in cells
            if (phase == "all" or r["phase"] == phase)
            and (r["arm"], r["score"]) == (arm, score)
        ]
        scores = {(r["receiver"], r["seed"], r["ligand"]): r["energy"] for r in sel}
        for rec in ("all", *RECEIVERS):
            aligned = [
                r
                for r in sel
                if (r["receiver"] == rec or rec == "all") and r["rmsd"] is not None
            ]
            blocks = {(r, s) for r, s, lig in scores if rec == "all" or rec == r}
            rank = []
            for a, b in itertools.combinations(sorted(series), 2):
                a, b = sorted(
                    (a, b), key=lambda l: float(series[l]["Kd_uM_derived_at_298.15K"])
                )
                if float(series[b]["Kd_uM_derived_at_298.15K"]) < 3 * float(
                    series[a]["Kd_uM_derived_at_298.15K"]
                ):
                    continue
                gap = median(scores[r, s, b] - scores[r, s, a] for r, s in blocks)
                rank.append(
                    {
                        "phase": phase,
                        "arm": arm,
                        "score": score,
                        "receiver": rec,
                        "strong": a,
                        "weak": b,
                        "predicted_gap": gap,
                        "correct": 1 if gap > 1e-9 else 0 if gap < -1e-9 else 0.5,
                    }
                )
            pairs.extend(rank)
            same_charge = [v for v in rank if "4" not in (v["strong"], v["weak"])]
            cross = [r for r in aligned if r["kind"] == "cross"]
            summary.append(
                {
                    "phase": phase,
                    "arm": arm,
                    "score": score,
                    "receiver": rec,
                    "pose_n": len(aligned),
                    "pose_hit": sum(r["rmsd"] <= 2 for r in aligned),
                    "oracle": sum(r["best_rmsd"] <= 2 for r in aligned),
                    "cross_n": len(cross),
                    "cross_hit": sum(r["rmsd"] <= 2 for r in cross),
                    "cross_oracle": sum(r["best_rmsd"] <= 2 for r in cross),
                    "pair_n": len(rank),
                    "pair_correct": sum(v["correct"] for v in rank),
                    "same_charge_n": len(same_charge),
                    "same_charge_correct": sum(v["correct"] for v in same_charge),
                }
            )
    for phase in ("first", "confirmation", "all"):
        pe = [
            r
            for r in pairs
            if r["phase"] == phase and r["receiver"] == "all" and r["score"] == "native"
        ]
        pp = {(r["strong"], r["weak"], r["arm"]): r["correct"] for r in pe}
        for dropped in ("none", *series):
            edges = [(a, b) for a, b, k in pp if k == "hard" and dropped not in (a, b)]
            hard = sum(pp[a, b, "hard"] for a, b in edges)
            soft = sum(pp[a, b, "soft"] for a, b in edges)
            omit.append(
                {
                    "phase": phase,
                    "omitted": dropped,
                    "pairs": len(edges),
                    "hard": hard,
                    "soft": soft,
                    "delta_percentage_points": 100 * (soft - hard) / len(edges),
                }
            )
    for phase in ("first", "confirmation"):
        sel = [
            r
            for r in cells
            if r["phase"] == phase and r["kind"] == "cross" and r["score"] == "native"
        ]
        by = {
            (r["receiver"], r["ligand"], r["seed"], r["arm"]): r["rmsd"] <= 2
            for r in sel
        }
        changes = [
            (v, by[rec, lig, s, "soft"])
            for (rec, lig, s, arm), v in by.items()
            if arm == "hard"
        ]
        paired.append(
            {
                "phase": phase,
                "both": sum(a and b for a, b in changes),
                "neither": sum(not a and not b for a, b in changes),
                "gain": sum(not a and b for a, b in changes),
                "loss": sum(a and not b for a, b in changes),
            }
        )
    for name, table in [
        ("soft_cells", cells),
        ("soft_pairs", pairs),
        ("soft_summary", summary),
        ("soft_leave_ligand_out", omit),
        ("soft_pose_switches", paired),
    ]:
        dump(out / (name + ".csv"), table)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("hard_results", type=Path)
    ap.add_argument("soft_results", type=Path)
    ap.add_argument("confirm_results", type=Path)
    ap.add_argument("out", type=Path)
    a = ap.parse_args()
    a.out.mkdir(exist_ok=True, parents=True)
    series = {r["id"]: r for r in read_csv("series.csv")}
    cells = [
        v
        for phase in SEEDS
        for v in evaluate(
            phase, a.hard_results, a.soft_results, a.confirm_results, a.out, series
        )
    ]
    summarize(cells, series, a.out)
