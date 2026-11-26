"""Test externally inferred water occupancy on published poses and same-assay ligands."""

import argparse
import csv
import hashlib
import itertools
import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from water_occupancy import atom_arrays, box, build_sites, read_models, score_components

DONORS = {
    "brd4": ("2OSS", "3ZYU", "4J0S", "4MR4"),
    "ca2": ("1CA2", "2CBA", "2ILI", "3KKX"),
}
RECS = {"brd4": ("3MXF", "4OGI", "4OGJ"), "ca2": ("6WQ4", "6WQ8")}
TARGET_OF = {"brd4": "brd4", "ca2": "ca2", "ca2_both": "ca2"}


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def brd_pools(brd):
    kd = {r["ligand"]: float(r["kd_uM"]) for r in read(brd / "cells.csv")}
    for (rec, lig, seed, arm), rows in itertools.groupby(
        sorted(read(brd / "modes.csv"), key=brd_key), key=brd_key
    ):
        yield (
            "brd4",
            rec,
            lig,
            seed,
            arm,
            kd[lig],
            [
                (
                    str(r["model"]),
                    brd / f"{rec}_{lig}_{seed}_{arm}.pdbqt",
                    int(r["model"]),
                    float(r["score"]),
                    float(r["rmsd"]),
                    None,
                )
                for r in rows
            ],
        )


def iterate(args):
    yield from brd_pools(args.brd_pools)
    ca = args.ca_pools
    ki = {r["lig"]: float(r["ki_nM"]) for r in read(ca / "cells.csv")}
    for (rec, lig, seed, arm), rows in itertools.groupby(
        sorted(read(ca / "modes.csv"), key=ca_key), key=ca_key
    ):
        yield (
            "ca2",
            rec,
            lig,
            seed,
            arm,
            ki[lig],
            [
                (
                    str(r["model"]),
                    ca / f"{rec}_{lig}_{seed}_{arm}.pdbqt",
                    int(r["model"]),
                    float(r["native"]),
                    float(r["rmsd"]),
                    None,
                )
                for r in rows
            ],
        )
    state = args.state_pools
    k = lambda r: (r["receptor"], r["ligand"], r["seed"])
    all_rows = [
        r
        for r in read(state / "candidates.csv")
        if r["slot"] in ("common", "second_reflected")
    ]
    for (rec, lig, seed), rows in itertools.groupby(sorted(all_rows, key=k), key=k):
        rows = list(rows)
        best = min(
            rows,
            key=lambda r: (float(r["total"]), r["slot"] != "common", int(r["model"])),
        )
        unbound = float(best["intra"])
        yield (
            "ca2_both",
            rec,
            lig,
            seed,
            "scaffold",
            ki[lig],
            [
                (
                    r["slot"] + ":" + r["model"],
                    state / f"{rec}_{lig}_{seed}_{r['slot']}.pdbqt",
                    int(r["model"]),
                    (float(r["total"]) - unbound) * float(r["flex_factor"]),
                    float(r["rmsd"]),
                    bool(int(r["same_sign_as_bound"])),
                )
                for r in rows
            ],
        )


def brd_key(r):
    return r["receptor"], r["ligand"], r["seed"], r["arm"]


def ca_key(r):
    return r["rec"], r["lig"], r["seed"], r["arm"]


@lru_cache(maxsize=16)
def parse_pool(path):
    return {
        i: atom_arrays(lines)[:2]
        for i, lines in read_models(path.read_text().splitlines()).items()
    }


def calculate(args):
    sites = {}
    donor_hash = {
        f.stem: hashlib.sha256(f.read_bytes()).hexdigest()
        for f in args.donors.glob("*.pdb")
        if any(f.stem in v for v in DONORS.values())
    }
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    (out / "donors.json").write_text(
        json.dumps(donor_hash, sort_keys=True, indent=2) + "\n"
    )
    for target, recs in RECS.items():
        work = args.brd if target == "brd4" else args.ca
        names = DONORS[target]
        for rec in recs:
            for leave in ("all", *names):
                donors = [
                    (d, "A", args.donors / (d + ".pdb")) for d in names if d != leave
                ]
                chain = "B" if rec == "4OGJ" else "A"
                prep = work / "prepared"
                evidence = build_sites(
                    work / "raw" / (rec + ".pdb"),
                    rec,
                    chain,
                    prep / (rec + "_receptor.pdbqt"),
                    *box(prep / (rec + "_box.txt")),
                    donors,
                )
                sites[target, rec, leave] = evidence["sites"]
                if leave == "all":
                    (out / (rec + "_sites.json")).write_text(
                        json.dumps(evidence, indent=2) + "\n"
                    )
    poses = []
    cells = []
    pairs = []
    summary = []
    for name, rec, lig, seed, arm, k, models in iterate(args):
        target = TARGET_OF[name]
        base = {"panel": name, "rec": rec, "lig": lig, "seed": int(seed), "arm": arm}
        for leave in ("all", *DONORS[target]):
            cand = []
            for mid, path, index, native, rms, sign in models:
                xyz, polar = parse_pool(path)[index]
                vals = score_components(sites[target, rec, leave], xyz, polar)
                row = {
                    **base,
                    "candidate": mid,
                    "label": k,
                    "rmsd": rms,
                    "sign": sign,
                    "native": native,
                    **{m: native + d for m, d in vals.items()},
                    **{"delta_" + m: d for m, d in vals.items()},
                }
                cand.append(row)
                if leave == "all":
                    poses.append(row)
            for method in (
                (
                    "native",
                    "exclusion",
                    "contact",
                    "joint",
                    "replace_any",
                    "replace_polar",
                )
                if leave == "all"
                else ("joint", "replace_polar")
            ):
                win = min(cand, key=lambda r: r[method])
                scores = [r for r in cand if r["rmsd"] <= 2 and r["sign"] is not False]
                false = [r for r in cand if r not in scores]
                separating = (
                    (min(r[method] for r in false) - min(r[method] for r in scores))
                    if scores and false
                    else None
                )
                cells.append(
                    {
                        **base,
                        "leave_out": leave,
                        "method": method,
                        "candidate": win["candidate"],
                        "rmsd": win["rmsd"],
                        "signed_hit": int(win in scores),
                        "oracle": int(bool(scores)),
                        "score": win[method],
                        "native_at_selected": win["native"],
                        "label": k,
                        "positive_gap": separating,
                    }
                )

    def key(r):
        return r["panel"], r["arm"], r["leave_out"], r["method"]

    for (panel, arm, leave, method), group in itertools.groupby(
        sorted(cells, key=key), key
    ):
        group = list(group)
        base = {"panel": panel, "arm": arm, "leave_out": leave, "method": method}
        edges = []
        for (rec, seed), block in itertools.groupby(
            sorted(group, key=lambda r: (r["rec"], r["seed"])),
            lambda r: (r["rec"], r["seed"]),
        ):
            for a, b in itertools.combinations(
                sorted(block, key=lambda r: r["label"]), 2
            ):
                if b["label"] < 3 * a["label"]:
                    continue
                delta = b["score"] - a["score"]
                native_delta = b["native_at_selected"] - a["native_at_selected"]
                edges.append(
                    {
                        **base,
                        "rec": rec,
                        "seed": seed,
                        "strong": a["lig"],
                        "weak": b["lig"],
                        "score_gap": delta,
                        "native_at_selected_gap": native_delta,
                        "credit": (1 if delta > 0 else 0.5 if delta == 0 else 0),
                    }
                )
        pairs.extend(edges)

        def pairkey(r):
            return r["strong"], r["weak"]

        gaps = [
            np.median([r["score_gap"] for r in g])
            for _, g in itertools.groupby(sorted(edges, key=pairkey), pairkey)
        ]
        summary.append(
            {
                **base,
                "cells": len(group),
                "top_2A": sum(r["rmsd"] <= 2 for r in group),
                "joint_top": sum(r["signed_hit"] for r in group),
                "joint_pool": sum(r["oracle"] for r in group),
                "median_rmsd": np.median([r["rmsd"] for r in group]),
                "pairs": len(edges),
                "credit": sum(r["credit"] for r in edges),
                "distinct_pairs": len(gaps),
                "distinct_credit": sum(
                    1 if v > 0 else 0.5 if v == 0 else 0 for v in gaps
                ),
            }
        )
    for file, data in (
        ("poses", poses),
        ("cells", cells),
        ("pairs", pairs),
        ("summary", summary),
    ):
        write(out / (file + ".csv"), data)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("donors", "brd", "ca", "brd_pools", "ca_pools", "state_pools", "out"):
        p.add_argument("--" + key.replace("_", "-"), type=Path, required=True)
    calculate(p.parse_args())
