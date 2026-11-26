"""Select poses by blind component scores, then compare paired RMSD and observed affinity."""

import argparse
import functools
import itertools
import operator
import statistics
from collections import defaultdict
from pathlib import Path

from run import read, write

METHODS = (
    "native",
    "coulomb",
    "polarization",
    "cycle",
    "replace_hbond",
    "cycle_raw_scale",
    "eps2",
    "eps8",
    "cycle_only",
    "cycle_at_native",
    "native_at_cycle",
    "retain_scaffold",
    "retain_scaffold_cycle",
)


def energy(r, zero, method):
    if method == "cycle_at_native":
        method = "cycle"
    if method == "native_at_cycle":
        method = "native"
    factor = float(r["factor"])
    native = factor * (float(r["inter"]) + float(r["intra"]) - zero)
    c, gb = float(r["coulomb"]), float(r["delta_gb"])
    correction = {
        "native": 0,
        "coulomb": c,
        "polarization": gb,
        "cycle": c + gb,
        "replace_hbond": c + gb - factor * float(r["hbond"]),
        "cycle_raw_scale": factor * (c + gb),
    }
    if method == "cycle_only":
        return c + gb
    if method in ("retain_scaffold", "retain_scaffold_cycle"):
        return (
            native
            + factor * float(r["prior"])
            + (c + gb if method.endswith("_cycle") else 0)
        )
    if method in ("eps2", "eps8"):
        eps = int(method[3:])
        return native + c * 4 / eps + gb * (1 / eps - 1 / 78.5) / (1 / 4 - 1 / 78.5)
    return native + correction[method]


def select(rows, prior):
    bycell = defaultdict(list)
    for r in rows:
        key = r["target"], r["receptor"], r["ligand"], r["seed"]
        if r["target"] == "hcaii":
            r["prior"] = prior[(*key[1:], r["slot"], r["model"])]
            for arm, second in (
                ("both_states", "second_reflected"),
                ("original_only", "second_original"),
            ):
                if r["slot"] in ("common", second):
                    bycell[(*key, arm)].append(r)
        else:
            bycell[(*key, r["slot"])].append(r)
    cells = []
    for (target, rec, lig, seed, arm), pool in sorted(bycell.items()):
        tie = lambda r: (r["slot"] != "common", int(r["model"]))
        native = min(
            pool, key=lambda r: (float(r["inter"]) + float(r["intra"]), *tie(r))
        )
        zero = float(native["intra"])
        hit = lambda r: float(r["rmsd"]) <= 2 and int(r["sign"])
        oracle = any(hit(r) for r in pool)
        for method in METHODS:
            if method.startswith("retain_scaffold") and target != "hcaii":
                continue
            chooser = {"cycle_at_native": "native", "native_at_cycle": "cycle"}.get(
                method, method
            )
            chosen = min(pool, key=lambda r: (energy(r, zero, chooser), *tie(r)))
            cells.append(
                {
                    "target": target,
                    "receptor": rec,
                    "ligand": lig,
                    "seed": seed,
                    "pool": arm,
                    "method": method,
                    "slot": chosen["slot"],
                    "model": chosen["model"],
                    "rmsd": chosen["rmsd"],
                    "sign": chosen["sign"],
                    "hit": int(hit(chosen)),
                    "oracle": int(oracle),
                    "exp": chosen["exp"],
                    "score": energy(chosen, zero, method),
                    "native_at_selected": energy(chosen, zero, "native"),
                }
            )
    return cells


def evaluate(cells):
    paired = defaultdict(list)
    grouped = defaultdict(list)
    for r in cells:
        key = r["target"], r["pool"], r["method"]
        paired[(*key, r["receptor"], r["seed"])].append(r)
        grouped[key].append(r)
    pairs = []
    pair_groups = defaultdict(list)
    for (target, arm, method, rec, seed), rows in sorted(paired.items()):
        for a, b in itertools.combinations(
            sorted((r for r in rows if r["exp"] != ""), key=lambda r: float(r["exp"])),
            2,
        ):
            if float(b["exp"]) / float(a["exp"]) < 3:
                continue
            diff = b["score"] - a["score"]
            q = {
                "target": target,
                "pool": arm,
                "method": method,
                "receptor": rec,
                "seed": seed,
                "strong": a["ligand"],
                "weak": b["ligand"],
                "diff": diff,
                "correct": 1.0 if diff > 0 else 0.5 if diff == 0 else 0.0,
            }
            pairs.append(q)
            pair_groups[target, arm, method].append(q)
    summary = []
    for (target, arm, method), rows in sorted(grouped.items()):
        p = pair_groups[target, arm, method]
        bypair = defaultdict(list)
        for q in p:
            bypair[q["strong"], q["weak"]].append(q["diff"])
        med = [statistics.median(d) for d in bypair.values()]
        summary.append(
            {
                "target": target,
                "pool": arm,
                "method": method,
                "cells": len(rows),
                "rmsd_top": sum(float(r["rmsd"]) <= 2 for r in rows),
                "strict_top": sum(r["hit"] for r in rows),
                "strict_oracle": sum(r["oracle"] for r in rows),
                "median_rmsd": statistics.median(float(r["rmsd"]) for r in rows),
                "pairs_correct": sum(q["correct"] for q in p),
                "pairs": len(p),
                "distinct_correct": sum(int(x > 0) + 0.5 * (x == 0) for x in med),
                "distinct": len(med),
            }
        )
    return pairs, summary


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results", type=Path)
    p.add_argument(
        "scaffold_prior", type=Path, help="Original ammonium rank_energies.csv"
    )
    a = p.parse_args()
    pr = {
        (r["receptor"], r["ligand"], r["seed"], r["slot"], r["model"]): float(
            r["prior"]
        )
        for r in read(a.scaffold_prior)
    }
    rows = functools.reduce(operator.iadd, [
            read(a.results / (k + ".csv"))
            for k in ("hcaii", "brd4", "t4l", "ca_holdout36", "ca_holdout20")
        ], [])
    cells = select(rows, pr)
    pairs, summary = evaluate(cells)
    for name, data in (("cells", cells), ("pairs", pairs), ("summary", summary)):
        write(a.results / (name + ".csv"), data)
