"""Recompute disjoint seed panels, substitution cycles, ligand-group uncertainty and preparation effects."""

import argparse
import math
from pathlib import Path
from statistics import median
import numpy as np

from experiment import ARMS, ROOT, SEEDS, table, write
from evaluate import RT, compare, pose_models
from repair import LIG

EXTRA = (161803, 141421)
SQUARE = {
    "Ai": "CHEMBL177142",
    "An": "CHEMBL176813",
    "Bi": "CHEMBL176744",
    "Bn": "CHEMBL369161",
}


def gap(cells, target, arm, method, ligs, coeff):
    by = {
        (r["receptor"], r["seed"], r["ligand"]): float(r["score"])
        for r in cells
        if (r["target"], r["arm"], r["method"]) == (target, arm, method)
    }
    cond = sorted(set((r, s) for r, s, k in by))
    values = [
        sum(c * by[r, s, k] for k, c in zip(ligs, coeff, strict=True)) for r, s in cond
    ]
    return median(values), min(values), max(values)


def unique_resampling(edges, ligs):
    ids = {s: i for i, s in enumerate(ligs)}
    idx = np.array([(ids[a], ids[b]) for a, b, d in edges], dtype=int)
    deltas = np.array([d for a, b, d in edges])
    # Vertex/bootstrap resamples entire chemical identities, preserving dependence of edges.
    rng = np.random.default_rng(226010)
    counts = rng.multinomial(len(ligs), np.full(len(ligs), 1 / len(ligs)), size=3000)
    weights = counts[:, idx[:, 0]] * counts[:, idx[:, 1]]
    total = weights.sum(axis=1)
    frac = (weights @ deltas)[total > 0] / total[total > 0]
    loo = [
        sum(d for a, b, d in edges if s not in (a, b))
        / sum(s not in (a, b) for a, b, d in edges)
        for s in ligs
        if any(s not in (a, b) for a, b, d in edges)
    ]
    return tuple(np.percentile(frac, [2.5, 97.5])), (min(loo), max(loo))


def evaluate(work, results, out):
    cells = table(results / "cells.csv")
    for r in cells:
        r["score"] = float(r["score"])
    metadata = [
        r
        for r in table(ROOT / "benchmarks/affinity_partition/fxa_series.csv")
        if int(r["included"])
    ]
    measures = {r["ligand"]: float(r["ki_nm"]) for r in metadata}
    labels = {
        "fxa": measures,
        "trypsin": {
            r["ligand"]: float(r["kd_uM"]) for r in table(work / "trypsin/series.csv")
        },
    }
    measuresize = {r["ligand"]: (r["heavy"], r["charge"]) for r in metadata}
    excluded = "CONH2"
    stats, cycles, pose_sum = [], [], []
    substitutions = (
        ("branch_A", ["An", "Ai"], [1, -1]),
        ("branch_B", ["Bn", "Bi"], [1, -1]),
        ("position_n", ["Bn", "An"], [1, -1]),
        ("position_i", ["Bi", "Ai"], [1, -1]),
        ("coupling", ["An", "Ai", "Bn", "Bi"], [1, -1, -1, 1]),
        ("diagonal", ["Bn", "Ai"], [1, -1]),
    )
    repaired = [
        dict(r) for r in cells if r["method"] == "release" and r["target"] == "fxa"
    ]
    for r in repaired:
        if r["ligand"] == LIG:
            path = (
                work / "repair" / f"{r['receptor']}_{LIG}_{r['seed']}_{r['arm']}.pdbqt"
            )
            energy, _ = next(pose_models(path))
            r["score"] = energy
    for split, selected in (
        ("initial", SEEDS),
        ("new", EXTRA),
        ("all", (*SEEDS, *EXTRA)),
    ):
        for prep, raw in (("original", cells), ("repaired", repaired)):
            selected_cells = [r for r in raw if int(r["seed"]) in selected]
            for rec in (None, "1FJS", "2P16"):
                own = [r for r in selected_cells if not rec or r["receptor"] == rec]
                pairs = compare(own, labels)
                for target in ("fxa", "trypsin"):
                    ligands = sorted(labels[target])
                    if not any(r["target"] == target for r in own):
                        continue
                    for method in ("release", "competition", "pair_rerank"):
                        if prep == "repaired" and method != "release":
                            continue
                        for subset in (
                            "all",
                            "matched",
                            "unambiguous",
                            "no_distortion",
                        ):
                            if (subset == "unambiguous" and target != "trypsin") or (
                                subset in ("matched", "no_distortion")
                                and target != "fxa"
                            ):
                                continue
                            keep = lambda r: (
                                r["target"] == target
                                and r["ratio"] >= 3
                                and r["method"] == method
                                and (
                                    subset != "matched"
                                    or measuresize[r["strong"]]
                                    == measuresize[r["weak"]]
                                )
                                and (
                                    subset != "unambiguous"
                                    or excluded not in (r["strong"], r["weak"])
                                )
                                and (
                                    subset != "no_distortion"
                                    or LIG not in (r["strong"], r["weak"])
                                )
                            )
                            base = {
                                (r["strong"], r["weak"]): r["credit"]
                                for r in pairs
                                if keep(r) and r["arm"] == "native"
                            }
                            if not base:
                                continue
                            for arm in ARMS:
                                edges = [
                                    r for r in pairs if keep(r) and r["arm"] == arm
                                ]
                                changes = [
                                    (
                                        r["strong"],
                                        r["weak"],
                                        r["credit"] - base[r["strong"], r["weak"]],
                                    )
                                    for r in edges
                                ]
                                shifts = [d for a, b, d in changes]
                                dat = dict(
                                    split=split,
                                    prep=prep,
                                    target=target,
                                    receptor=rec or "both",
                                    method=method,
                                    subset=subset,
                                    arm=arm,
                                    correct=sum(r["credit"] for r in edges),
                                    pairs=len(edges),
                                    rescues=sum(x > 0 for x in shifts),
                                    losses=sum(x < 0 for x in shifts),
                                )
                                if (
                                    subset == "all"
                                    and rec is None
                                    and prep == "original"
                                    and method == "release"
                                ):
                                    ci, loo = unique_resampling(changes, ligands)
                                    dat.update(
                                        boot_low=ci[0],
                                        boot_high=ci[1],
                                        loo_low=loo[0],
                                        loo_high=loo[1],
                                    )
                                else:
                                    dat.update(
                                        boot_low="",
                                        boot_high="",
                                        loo_low="",
                                        loo_high="",
                                    )
                                stats.append(dat)
                    if rec or prep != "original":
                        continue
                    for arm in ARMS:
                        rows = [
                            r
                            for r in own
                            if (r["target"], r["arm"], r["method"])
                            == (target, arm, "release")
                        ]
                        bens = [r for r in rows if r["heavy_rms"] != ""]
                        pose_sum.append(
                            dict(
                                split=split,
                                target=target,
                                arm=arm,
                                cells=len(rows),
                                top_salt=sum(int(r["bidentate"]) for r in rows),
                                pool_salt=sum(int(r["oracle_bidentate"]) for r in rows),
                                crystal_cells=len(bens),
                                top2A=sum(float(r["heavy_rms"]) <= 2 for r in bens),
                                pool2A=sum(float(r["oracle_rms"]) <= 2 for r in bens),
                            )
                        )
            if prep == "original":
                for name, ids, coeff in substitutions:
                    ligs = [SQUARE[k] for k in ids]
                    expected = RT * sum(
                        c * math.log(measures[l]) for l, c in zip(ligs, coeff)
                    )
                    for arm in ARMS:
                        mid, lo, hi = gap(
                            selected_cells, "fxa", arm, "release", ligs, coeff
                        )
                        cycles.append(
                            dict(
                                split=split,
                                effect=name,
                                arm=arm,
                                observed_kcal=expected,
                                median_kcal=mid,
                                min_kcal=lo,
                                max_kcal=hi,
                            )
                        )
    write(out / "contrast.csv", stats)
    write(out / "substitution_cycle.csv", cycles)
    write(out / "pose_summary.csv", pose_sum)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "results", "out"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    evaluate(a.work, a.results, a.out)
