"""Continuous ligand pose errors from the frozen 7ca6de9 studies.

The unit is one returned pool for a specified query, receptor, budget and seed.
Use --ecdf-out for exact step curves, with no binning or upper cutoff.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
OH = BENCH / "ligand_hydroxyl"
KEY = ["study", "task", "budget", "cap", "target", "query", "receiver", "seed"]
ARM = ["pool", "method"]
ERRORS = ["top1_A", "best5_A", "oracle_A", "regret_A"]
INPUT = KEY + ARM + ["reference_pool", "reference_method", "n_poses", *ERRORS[:3]]


def score_pools(frame, methods, default_pool):
    """One record per run and score. Native ties keep released pose order."""
    output = []
    for key, poses in frame.groupby(KEY + ["pool"], sort=False):
        poses = poses.sort_values("rank")
        if poses["rank"].tolist() != list(range(1, len(poses) + 1)):
            raise ValueError(f"Discontinuous output mode indices in {key}")
        for method, score_column in methods.items():
            if key[-1] == "upstream" and method != "native":
                continue
            if not np.isfinite(poses[[score_column, "rmsd"]]).all().all():
                raise ValueError(f"Missing score/RMSD in {key}, {method}")
            ordered = (
                poses
                if method == "native"
                else poses.sort_values([score_column, "rank"])
            )
            errors = ordered["rmsd"]
            output.append(
                dict(
                    zip(KEY + ["pool"], key),
                    method=method,
                    n_poses=len(poses),
                    top1_A=errors.iloc[0],
                    best5_A=errors.iloc[:5].min(),
                    oracle_A=errors.min(),
                    reference_pool=default_pool,
                    reference_method="native",
                )
            )
    return pd.DataFrame(output)


def hydroxyl_runs():
    archive = pd.read_csv(OH / "results/poses.csv")
    archive["study"] = "archive"
    archive["task"] = archive["panel"]
    archive["budget"] = archive["design"]
    archive["pool"] = "archived"
    methods = {"native": "native", "bonded": "bonded"}
    yield score_pools(archive, methods, "archived")

    transfer = pd.read_csv(OH / "transfer/results/poses.csv")
    transfer = transfer.rename(
        columns={"model": "rank", "engine": "pool", "kind": "task"}
    )
    transfer["study"] = "transfer"
    transfer["budget"] = transfer["cap"].astype(str)
    yield score_pools(transfer, methods, "pro")

    sugar = pd.read_csv(OH / "glucose/results/poses.csv").rename(
        columns={"engine": "pool"}
    )
    sugar["study"] = "glucose"
    sugar["task"] = "self"
    sugar["budget"] = sugar["cap"].astype(str)
    sugar[["receiver", "query"]] = sugar["case"].str.split("-", expand=True)
    sugar["target"] = sugar["receiver"].map(
        {"1V4S": "GCK", "3IDH": "GCK", "4ZW9": "GLUT3"}
    )
    yield score_pools(
        sugar,
        {
            "native": "native",
            "bonded": "donor",
            "bidirectional": "bidirectional",
            "lone_pairs": "lone_pairs",
        },
        "pro",
    )


def macrocycle_runs():
    root = BENCH / "macrocycle_pose_distance"
    x = pd.concat(
        [pd.read_csv(root / run / "structures.csv") for run in ("standard", "higher")]
    )
    x = x.rename(columns={"arm": "pool", "modes": "n_poses"})
    return x.assign(
        study="macrocycle",
        task="self",
        budget=x.cap.astype(str),
        target=x.pdb.map({"4EGK": "HSP90a", "3VHA": "HSP90a", "2J9M": "CDK2"}),
        receiver=x.pdb,
        query=x.pdb,
        method="native",
        reference_pool="parent",
        reference_method="native",
    )


def lipid_runs():
    root = OH / "glucose/lipid/results"
    x = pd.read_csv(root / "cells.csv").rename(
        columns={
            "score": "method",
            "code": "query",
            "modes": "n_poses",
            "top1": "top1_A",
            "top5": "best5_A",
            "oracle": "oracle_A",
        }
    )
    near = pd.read_csv(root / "distant_control.csv").rename(
        columns={
            "modes": "n_poses",
            "top1": "top1_A",
            "top5": "best5_A",
            "oracle": "oracle_A",
        }
    )
    near["query"] = near["case"].str.split("-").str[1]
    x = pd.concat([x, near.assign(pool="near_only", method="reported")])
    return x.assign(
        study="lipid",
        task=x["query"],
        budget=x.cap.astype(str),
        target="GLUT3",
        receiver="4ZW9",
        reference_pool="dry",
        reference_method=np.where(x.method == "reported", "reported", "dry"),
    )


def pair_runs(runs=None):
    if runs is None:
        runs = pd.concat([*hydroxyl_runs(), macrocycle_runs(), lipid_runs()])
    missing = set(INPUT) - set(runs)
    if missing:
        raise ValueError(f"Missing run columns: {sorted(missing)}")
    runs = runs[INPUT].copy().reset_index(drop=True)
    if runs.empty or runs.isna().any().any():
        raise ValueError("Expected non-empty run records with all required values")
    numeric = ["cap", "seed", "n_poses", *ERRORS[:3]]
    labels = [c for c in INPUT if c not in numeric]
    for c in labels:
        runs[c] = runs[c].astype(str).str.strip()
        if runs[c].eq("").any():
            raise ValueError(f"Empty label in {c}")
    for c in numeric:
        runs[c] = pd.to_numeric(runs[c], errors="raise")
    if not np.isfinite(runs[numeric]).all().all():
        raise ValueError("Nonfinite RMSD, seed, cap, or pool size")
    for c in ["n_poses", "cap", "seed"]:
        if (runs[c] % 1 != 0).any():
            raise ValueError(f"Expected integer {c}")
    if runs.duplicated(KEY + ARM).any():
        raise ValueError("Duplicate run and method")
    if (
        (runs.top1_A < runs.best5_A)
        | (runs.best5_A < runs.oracle_A)
        | (runs.oracle_A < 0)
        | (runs.n_poses < 1)
        | (runs.cap < 0)
        | ((runs.n_poses <= 5) & (runs.best5_A != runs.oracle_A))
        | ((runs.n_poses == 1) & (runs.top1_A != runs.oracle_A))
    ).any():
        raise ValueError("Invalid pool size, budget, or incompatible RMSDs")
    pool_counts = runs.groupby(KEY + ["pool"])["n_poses"].nunique()
    if (pool_counts != 1).any():
        raise ValueError(
            "Methods on the same pool must report the same number of poses"
        )
    runs["regret_A"] = runs.top1_A - runs.oracle_A
    ref = runs[KEY + ARM + ERRORS].rename(
        columns={
            "pool": "reference_pool",
            "method": "reference_method",
            **{e: "baseline_" + e for e in ERRORS},
        }
    )
    pairs = runs.merge(
        ref,
        on=KEY + ["reference_pool", "reference_method"],
        how="left",
        validate="many_to_one",
    )
    for e in ERRORS:
        pairs["delta_" + e] = pairs[e] - pairs.pop("baseline_" + e)
    if pairs.isna().any().any():
        raise ValueError("An arm lacks a matching run under its reference method")
    same_pool = pairs.pool == pairs.reference_pool
    if (pairs.loc[same_pool, "delta_oracle_A"] != 0).any():
        raise ValueError("Pool coverage changed during pure re-ranking")
    pairs["comparison"] = np.select(
        [~same_pool, pairs.method != pairs.reference_method],
        ["separate_search", "fixed_pool"],
        default="reference",
    )
    return pairs.sort_values(KEY + ARM).reset_index(drop=True)


def balanced_weights(frame):
    # Successively sample a target, a deposited query within target, a receptor
    # within query, then a repeat, uniformly at each level.
    heads = []
    w = np.ones(len(frame))
    for column in ("target", "query", "receiver", "seed"):
        n = (
            frame.groupby(heads)[column].transform("nunique")
            if heads
            else frame[column].nunique()
        )
        w /= np.asarray(n)
        heads.append(column)
    if not np.isclose(w.sum(), 1):
        raise ValueError("Grouping failed to account for every run exactly once")
    return w


def quantiles(x, w, probs):
    """Inverse weighted empirical CDF; average endpoints at exact mass boundaries (type 2)."""
    order = np.argsort(x, kind="stable")
    x = np.asarray(x)[order]
    cdf = np.cumsum(w[order])
    result = []
    for p in probs:
        j = int(np.searchsorted(cdf, p - 1.0e-12))
        if j + 1 < len(x) and abs(cdf[j] - p) < 1.0e-12:
            result.append((x[j] + x[j + 1]) / 2)
        else:
            result.append(x[j])
    return result


def summary_block(frame, label, target, w, ecdf):
    row = {
        "scope": label,
        "target": target,
        "n_targets": frame.target.nunique(),
        "n_queries": frame[["target", "query"]].drop_duplicates().shape[0],
        "n_query_receiver": frame[["target", "query", "receiver"]]
        .drop_duplicates()
        .shape[0],
        "n_runs": len(frame),
        "min_pool_size": frame.n_poses.min(),
        "max_pool_size": frame.n_poses.max(),
    }
    for e in ERRORS:
        x = frame[e].to_numpy()
        row["mean_" + e] = np.dot(w, x)
        q25, q50, q75, q90 = quantiles(x, w, (0.25, 0.5, 0.75, 0.9))
        row["median_" + e] = q50
        if e == "top1_A":
            row.update(q25_top1_A=q25, q75_top1_A=q75, q90_top1_A=q90)
        delta = frame["delta_" + e].to_numpy()
        row["mean_delta_" + e] = np.dot(w, delta)
        row["median_delta_" + e] = quantiles(delta, w, [0.5])[0]
        if e == "top1_A":
            row.update(
                improved_weight=w[delta < 0].sum(),
                tied_weight=w[delta == 0].sum(),
                worsened_weight=w[delta > 0].sum(),
            )
        if e != "regret_A":
            row["success2_" + e] = w[x <= 2].sum()
            if ecdf is not None:
                values, inverse = np.unique(x, return_inverse=True)
                cdf = np.cumsum(np.bincount(inverse, weights=w))
                ecdf.extend(
                    {
                        "scope": label,
                        "target": target,
                        "metric": e,
                        "threshold_A": v,
                        "cdf": f,
                    }
                    for v, f in zip(values, cdf)
                )
    coverage = row["success2_oracle_A"]
    row["conditional_top1_at_2A"] = (
        row["success2_top1_A"] / coverage if coverage else np.nan
    )
    return row


def summarize(runs, include_curves=False):
    summaries, curves = [], []

    def emit(frame, meta, label, target, w):
        new_curves = [] if include_curves else None
        summaries.append(
            dict(**meta, **summary_block(frame, label, target, w, new_curves))
        )
        if new_curves:
            curves.extend(dict(**meta, **point) for point in new_curves)

    group_cols = [
        "study",
        "task",
        "budget",
        *ARM,
        "comparison",
        "reference_pool",
        "reference_method",
    ]
    runs["weight_in_stratum"] = np.nan
    runs["weight_in_target"] = np.nan
    for key, frame in runs.groupby(group_cols):
        meta = dict(zip(group_cols, key))
        w = balanced_weights(frame)
        runs.loc[frame.index, "weight_in_stratum"] = w
        emit(frame, meta, "target_balanced", "ALL", w)
        for target, subset in frame.groupby("target"):
            w = balanced_weights(subset)
            runs.loc[subset.index, "weight_in_target"] = w
            emit(subset, meta, "within_target", target, w)
    # Matches the earlier 136-run headline only. Multiple protocols and repeated
    # queries share these runs; it is not a target-level generalization estimate.
    for method, frame in runs[runs.study == "archive"].groupby("method"):
        meta = {
            "study": "archive",
            "task": "all_archived_runs",
            "budget": "mixed",
            "pool": "archived",
            "method": method,
            "comparison": frame.comparison.iloc[0],
            "reference_pool": "archived",
            "reference_method": "native",
        }
        emit(
            frame,
            meta,
            "unweighted_descriptive",
            "ALL",
            np.full(len(frame), 1 / len(frame)),
        )
    return pd.DataFrame(summaries), pd.DataFrame(curves)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=HERE / "results")
    p.add_argument("--ecdf-out", type=Path, help="Optional complete weighted ECDF CSV")
    p.add_argument(
        "--runs",
        type=Path,
        help="Evaluate a canonical run-level CSV; required columns: "
        + ", ".join(INPUT),
    )
    a = p.parse_args()
    data = None
    if a.runs is not None:
        labels = set(INPUT) - {"cap", "seed", "n_poses", *ERRORS[:3]}
        data = pd.read_csv(
            a.runs,
            dtype={c: str for c in labels},
            keep_default_na=False,
            float_precision="round_trip",
        )
    runs = pair_runs(data)
    summary, curves = summarize(runs, a.ecdf_out is not None)
    a.output.mkdir(parents=True, exist_ok=True)
    runs.to_csv(a.output / "pose_runs.csv", index=False, float_format="%.17g")
    summary.to_csv(a.output / "pose_summary.csv", index=False, float_format="%.12g")
    if a.ecdf_out is not None:
        a.ecdf_out.parent.mkdir(parents=True, exist_ok=True)
        curves.to_csv(a.ecdf_out, index=False, float_format="%.12g")


if __name__ == "__main__":
    main()
