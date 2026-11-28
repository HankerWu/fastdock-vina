"""Produce per-stratum rank, pair and actual-cost ledgers from the raw engine receipts."""

import argparse
from collections import defaultdict
import csv
from decimal import Decimal
from itertools import combinations
import json
from math import log10
from pathlib import Path
from statistics import median
import sys

import numpy as np
from scipy.stats import spearmanr
from run import ENGINES, HERE, STAGES, cells

sys.path.insert(0, str(HERE.parents[1] / "prediction_metrics"))
from affinity import rank_correlations


def write(path, rows):
    if rows:
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=rows[0])
            w.writeheader()
            w.writerows(rows)


def describe(compound):
    lines = compound.read_text().splitlines()
    heavy = sum(
        l.split()[-1] not in {"H", "HD", "HS"} for l in lines if l.startswith(("ATOM ", "HETATM"))
    )
    (tor,) = [int(l.split()[1]) for l in lines if l.startswith("TORSDOF ")]
    return heavy, tor


def analyze(work, out):
    out.mkdir(exist_ok=True, parents=True)
    summary, pred, decisions, costs, runs = [], [], [], [], []
    common = work / "common_cost.json"
    trials = [(stage, engine, stage) for stage in STAGES for engine in ENGINES]
    if common.is_file():
        timing = json.loads(common.read_text())
        for item in timing["attempts"]:
            trials.append((item["stage"], item["engine"], item["stage"]))
        trials.extend(("common", e, source) for e, source in timing["final"].items())
    for stage, engine, storage_stage in trials:
        groups, samples = defaultdict(dict), []
        for state, row, seed in cells(stage):
            base = dict(
                stage=stage,
                engine=engine,
                storage_stage=storage_stage,
                target=row["target"],
                assay=row["assay"],
                endpoint=row["endpoint"],
                receptor=state,
                compound=row["compound"],
                seed=seed,
            )
            p = work / "outputs" / engine / storage_stage / f"{state}_{row['compound']}_{seed}.json"
            r = json.loads(p.read_text()) if p.exists() else {}
            selected = {
                k: r.get(k)
                for k in (
                    "exit_status",
                    "valid",
                    "top_score",
                    "cap",
                    "modes",
                    "wall_s",
                    "user_s",
                    "sys_s",
                    "maxrss_kib",
                    "hunt_evals",
                    "second_evals",
                    "explicit_evals",
                    "late_prescore",
                    "output_sha256",
                )
            }
            runs.append({**base, **selected})
            if r:
                samples.append(r)
            identity = (row["target"], row["assay"], row["endpoint"], state)
            if row["compound"] not in groups[identity]:
                heavy, torsions = describe(work / "prepared" / (row["compound"] + ".pdbqt"))
                groups[identity][row["compound"]] = dict(
                    compound=row["compound"],
                    value_nM=float(row["value_nM"]),
                    name=row["name"],
                    formal_charge=int(row["formal_charge"]),
                    heavy_atoms=heavy,
                    torsions=torsions,
                    scaffold=row["scaffold"],
                    values=[],
                )
            groups[identity][row["compound"]]["values"].append(
                r["top_score"] if r.get("valid") else None
            )
        cost = dict(
            stage=stage,
            engine=engine,
            expected=sum(1 for _ in cells(stage)),
            launches=len(samples),
            valid=sum(r["valid"] for r in samples),
            cpu_s=sum(r["user_s"] + r["sys_s"] for r in samples),
            sum_wall_s=sum(r["wall_s"] for r in samples),
            peak_rss_kib=max((r["maxrss_kib"] for r in samples), default=None),
        )
        for k in ("hunt_evals", "second_evals", "explicit_evals", "late_prescore"):
            cost[k] = (
                sum(r[k] for r in samples)
                if samples and all(r[k] is not None for r in samples)
                else None
            )
        costs.append(cost)
        for (target, assay, endpoint, state), compounds in sorted(groups.items()):
            key = dict(
                stage=stage, engine=engine, target=target, assay=assay, endpoint=endpoint, receptor=state
            )
            available = []
            for m in compounds.values():
                scores = m.pop("values")
                m["score"] = round(median(scores), 4) if all(s is not None for s in scores) else None
                pred.append(
                    {
                        **key,
                        **m,
                        "complete_seeds": sum(s is not None for s in scores),
                        "expected_seeds": len(scores),
                    }
                )
                if m["score"] is not None:
                    available.append(m)
            pK = np.array([-log10(m["value_nM"] * 1e-9) for m in available])
            strength = np.array([-m["score"] for m in available])
            rho, tau = rank_correlations(pK, strength)
            result = dict(
                **key,
                expected_compounds=len(compounds),
                observed_compounds=len(available),
                rho=rho,
                tau_b=tau,
            )
            for prop in ("heavy_atoms", "torsions", "formal_charge"):
                x = np.array([m[prop] for m in available])
                for tag, response in (("predicted", strength), ("measured", pK)):
                    result[f"{prop}_{tag}_rho"] = (
                        spearmanr(x, response).statistic
                        if len(x) > 1 and np.ptp(x) > 0 and np.ptp(response) > 0
                        else None
                    )
            edges = []
            for first, second in combinations(compounds.values(), 2):
                a, b = sorted((first, second), key=lambda m: m["value_nM"])
                if a["value_nM"] == b["value_nM"]:
                    continue
                fold = Decimal(str(b["value_nM"])) / Decimal(str(a["value_nM"]))
                delta = (
                    None
                    if a["score"] is None or b["score"] is None
                    else round(b["score"] - a["score"], 7)
                )
                credit = None if delta is None else 1 if delta > 0 else 0 if delta < 0 else 0.5
                edge = dict(
                    **key,
                    stronger=a["compound"],
                    weaker=b["compound"],
                    fold=fold,
                    weaker_minus_stronger=delta,
                    credit=credit,
                    same_scaffold=a["scaffold"] == b["scaffold"],
                )
                edges.append(edge)
            decisions.extend(edges)
            for threshold in (1, 3, 10):
                expected = [edge for edge in edges if edge["fold"] >= threshold]
                obs = [edge for edge in expected if edge["credit"] is not None]
                result[f"expected_{threshold}x"] = len(expected)
                result[f"pairs_{threshold}x"] = len(obs)
                result[f"correct_{threshold}x"] = sum(edge["credit"] for edge in obs)
            summary.append(result)
    for name, table in (
        ("assays", summary),
        ("ligands", pred),
        ("pairs", decisions),
        ("costs", costs),
        ("runs", runs),
    ):
        write(out / f"{name}.csv", table)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    analyze(args.work, args.out)
