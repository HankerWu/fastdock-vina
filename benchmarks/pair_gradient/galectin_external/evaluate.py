"""Evaluate native first-pose predictions against the independent same-assay series."""

import argparse
import csv
import itertools
import math
import re
from decimal import Decimal
from pathlib import Path

from panel import write
from scipy.stats import kendalltau, spearmanr

RECEIVERS = ("1KJL", "3T1L")
COUNTS = ("hunt_evals", "second_evals", "explicit_evals")


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def evaluate(work):
    compounds = {int(r["compound"]): r for r in read(work / "manifest.csv")}
    fixed = read(work / "fixed.csv")
    parent = [r for r in fixed if r["arm"] == "parent"]
    methods = {
        "parent": parent,
        "linear": [r for r in fixed if r["arm"] == "linear"],
        "transported": read(work / "transported.csv"),
    }
    if (work / "matched.csv").is_file():
        methods["matched"] = read(work / "matched.csv")
    scores = {}
    costs = []
    pairs = []
    summary = []
    influences = []
    for method, rows in methods.items():
        for r in rows:
            stem = work / "docking" / "raw" / r["key"]
            rawscore = float(
                re.search(
                    r"^REMARK VINA RESULT:\s*([-.\d]+)",
                    stem.with_suffix(".pdbqt").read_text(),
                    re.MULTILINE,
                ).group(1)
            )
            log = stem.with_suffix(".log").read_text()
            actual = dict(re.findall(r"(\w+)=(\d+)", log))
            assert rawscore == float(r["first_score"]) and all(
                int(actual[k]) == int(r[k]) for k in COUNTS
            )
            scores[method, r["receptor"], int(r["compound"])] = rawscore
        for rec in (*RECEIVERS, "both"):
            selected = [r for r in rows if rec == "both" or r["receptor"] == rec]
            costs.append(
                {
                    "method": method,
                    "receiver": rec,
                    "cells": len(selected),
                    "elapsed_s": sum(int(r["elapsed_ns"]) for r in selected) / 1e9,
                    "joint": sum(sum(int(r[k]) for k in COUNTS) for r in selected),
                    "local_calls": sum(
                        int(r["hunt_calls"])
                        + int(r["second_calls"])
                        + int(r["explicit_calls"])
                        for r in selected
                    ),
                }
            )
        for i in compounds:
            scores[method, "consensus", i] = round(
                sum(scores[method, s, i] for s in RECEIVERS) / len(RECEIVERS), 4
            )
    for rec in (*RECEIVERS, "consensus"):
        for i, r in compounds.items():
            scores["size", rec, i] = -int(r["heavy"])
        for method in [*methods, "size"]:
            for i, j in itertools.combinations(compounds, 2):
                ri, rj = compounds[i], compounds[j]
                ki, kj = float(ri["kd_uM"]), float(rj["kd_uM"])
                di, dj = Decimal(ri["kd_uM"]), Decimal(rj["kd_uM"])
                if max(di, dj) < 3 * min(di, dj):
                    continue
                strong, weak = (i, j) if ki < kj else (j, i)
                gap = scores[method, rec, weak] - scores[method, rec, strong]
                pairs.append(
                    {
                        "receiver": rec,
                        "method": method,
                        "strong": strong,
                        "weak": weak,
                        "kd_ratio": max(ki, kj) / min(ki, kj),
                        "observed_ddG_over_RT": abs(math.log(ki / kj)),
                        "predicted_score_gap": gap,
                        "credit": 1 if gap > 0 else 0.5 if gap == 0 else 0,
                        "matched_ha_charge": ri["heavy"] == rj["heavy"]
                        and ri["charge"] == rj["charge"],
                    }
                )
            group = [r for r in pairs if r["receiver"] == rec and r["method"] == method]
            for subset in ("all", "equal_heavy_charge", "without_bivalent"):
                ids = [
                    i
                    for i in compounds
                    if subset != "without_bivalent" or i not in {45, 46}
                ]
                selected = [
                    r
                    for r in group
                    if r["strong"] in ids
                    and r["weak"] in ids
                    and (subset != "equal_heavy_charge" or r["matched_ha_charge"])
                ]
                kd = [float(compounds[i]["kd_uM"]) for i in ids]
                vals = [scores[method, rec, i] for i in ids]
                rho, tau = spearmanr(vals, kd).statistic, kendalltau(vals, kd).statistic
                summary.append(
                    {
                        "receiver": rec,
                        "method": method,
                        "subset": subset,
                        "n_compounds": len(ids),
                        "n_pairs": len(selected),
                        "credit": sum(r["credit"] for r in selected),
                        "spearman": rho if subset != "equal_heavy_charge" else "",
                        "kendall_tau_b": tau if subset != "equal_heavy_charge" else "",
                    }
                )
    for state in (*RECEIVERS, "consensus"):
        for subset in ("all", "equal_heavy_charge"):
            by_method = {
                m: [
                    r
                    for r in pairs
                    if r["receiver"] == state
                    and r["method"] == m
                    and (subset == "all" or r["matched_ha_charge"])
                ]
                for m in methods
            }
            for m in sorted(set(methods) - {"parent"}):
                for dropped in compounds:

                    def credit(rows, omitted):
                        return sum(
                            x["credit"]
                            for x in rows
                            if omitted not in (x["strong"], x["weak"])
                        )

                    influences.append(
                        {
                            "method": m,
                            "receiver": state,
                            "subset": subset,
                            "omitted": dropped,
                            "delta_credit": credit(by_method[m], dropped)
                            - credit(by_method["parent"], dropped),
                        }
                    )
    write(
        work / "scores.csv",
        [
            {"method": m, "receiver": r, "compound": i, "score": v}
            for (m, r, i), v in scores.items()
        ],
    )
    for filename, rows in [
        ("costs", costs),
        ("pairs", pairs),
        ("summary", summary),
        ("influence", influences),
    ]:
        write(work / (filename + ".csv"), rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    evaluate(p.parse_args().work)
