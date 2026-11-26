"""Independently rescore retained sites per BRD4 pose and audit matched rank changes."""

import argparse
import json
import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path

from check import atom_roles, blocks, scalar
from evaluate import brd_pools, read, write


def run(out, pools):
    key = lambda r: (r["rec"], r["lig"], str(r["seed"]), r["arm"])
    sites = {
        rec: json.loads((out / f"{rec}_sites.json").read_text())["sites"]
        for rec in ("3MXF", "4OGI", "4OGJ")
    }
    witness = {
        (r["rec"], int(r["index"]))
        for r in read(out / "witness_sites.csv")
        if r["displacement_witnesses"]
    }
    oracle = {
        (r["rec"], r["query"], int(r["site_index"]))
        for r in read(out / "query_site_audit.csv")
        if r["crystal_displaces"] == "1"
    }
    saved = {
        "dry": {
            key(r): (r["candidate"], float(r["score"]))
            for r in read(out / "cells.csv")
            if r["panel"] == "brd4"
            and r["method"] == "native"
            and r["leave_out"] == "all"
        },
        "joint": {
            key(r): (r["candidate"], float(r["score"]))
            for r in read(out / "cells.csv")
            if r["panel"] == "brd4"
            and r["method"] == "joint"
            and r["leave_out"] == "all"
        },
        "witness": {
            key(r): (r["candidate"], float(r["score"]))
            for r in read(out / "witness_cells.csv")
        },
        "oracle": {
            key(r): (r["candidate"], float(r["oracle_score"]))
            for r in read(out / "oracle_cells.csv")
            if r["panel"] == "brd4"
        },
    }
    cases = []
    for _, rec, lig, seed, arm, k, models in brd_pools(pools):
        cell = {"rec": rec, "lig": lig, "seed": seed, "arm": arm, "Kd_uM": k}
        coords = blocks(models[0][1])
        candidates = {method: [] for method in saved}
        retained = {
            "joint": sites[rec],
            "witness": [s for i, s in enumerate(sites[rec]) if (rec, i) not in witness],
            "oracle": [
                s for i, s in enumerate(sites[rec]) if (rec, lig, i) not in oracle
            ],
        }
        for candidate, _, index, dry, rms, _ in models:
            xyz, roles = atom_roles(coords[index])
            candidates["dry"].append((dry, candidate, rms))
            for method, allowed in retained.items():
                corrected = dry + float(scalar(xyz, roles, allowed)[2])
                candidates[method].append((corrected, candidate, rms))
        for method, scored in candidates.items():
            score, model, rms = min(scored, key=lambda r: r[0])
            expected, refscore = saved[method][key(cell)]
            if model != expected or not math.isclose(
                score, refscore, abs_tol=1e-10, rel_tol=0
            ):
                raise AssertionError((cell, method, model, score, expected, refscore))
            cell.update(
                {
                    method + "_model": model,
                    method + "_score": score,
                    method + "_rmsd": rms,
                    method + "_hit": int(rms <= 2),
                }
            )
        cases.append(cell)
    write(out / "mask_comparison_cells.csv", cases)
    grouped = defaultdict(list)
    for cell in cases:
        grouped[cell["rec"], cell["seed"], cell["arm"]].append(cell)
    pairs = []
    for (rec, seed, arm), block in grouped.items():
        for strong, weak in combinations(sorted(block, key=lambda r: r["Kd_uM"]), 2):
            if weak["Kd_uM"] < 3 * strong["Kd_uM"]:
                continue
            row = {
                "rec": rec,
                "seed": seed,
                "arm": arm,
                "strong": strong["lig"],
                "weak": weak["lig"],
            }
            for method in saved:
                gap = weak[method + "_score"] - strong[method + "_score"]
                row[method + "_gap"] = gap
                row[method + "_credit"] = int(gap > 0) + 0.5 * (gap == 0)
            pairs.append(row)
    write(out / "mask_comparison_pairs.csv", pairs)
    print(
        f"Rescored retained site lists independently for {len(cases)} BRD4 cells and {len(pairs)} paired measurements."
    )
    for arm in ("native", "transfer"):
        c = [r for r in cases if r["arm"] == arm]
        p = [r for r in pairs if r["arm"] == arm]
        for m in saved:
            print(
                arm,
                m,
                "top",
                sum(r[m + "_hit"] for r in c),
                "pair credit",
                sum(r[m + "_credit"] for r in p),
            )
        for f, t in (("dry", "joint"), ("joint", "witness"), ("joint", "oracle")):
            print(
                arm,
                f,
                t,
                "top gains/losses",
                sum(r[t + "_hit"] > r[f + "_hit"] for r in c),
                sum(r[t + "_hit"] < r[f + "_hit"] for r in c),
                "pair gains/losses",
                sum(r[t + "_credit"] > r[f + "_credit"] for r in p),
                sum(r[t + "_credit"] < r[f + "_credit"] for r in p),
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--brd-pools", required=True, type=Path)
    a = parser.parse_args()
    run(a.out, a.brd_pools)
