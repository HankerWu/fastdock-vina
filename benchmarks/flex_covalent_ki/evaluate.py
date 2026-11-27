"""Audit coverage, paired same-assay ranks and ligand-cluster bootstrap intervals."""

import argparse
import csv
import itertools
import math
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau, spearmanr

from run import SERIES, RECS, SEEDS, write


def ranking(y, x):
    # Average ranks handle ties; predicted strength increases with -Vina score.
    a, b = np.triu_indices(len(y), 1)
    dy = y[a] - y[b]
    dx = np.round(x[a] - x[b], 10)
    choose = np.abs(dy) >= np.log10(3) - 1e-12
    credits = 0.5 + 0.5 * np.sign(dy) * np.sign(dx)
    return np.asarray(
        [
            spearmanr(x, y).statistic,
            kendalltau(x, y, variant="b").statistic,
            credits[choose].mean() if choose.any() else np.nan,
        ]
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("runs", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--resamples", type=int, default=2000)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    runs = list(csv.DictReader(args.runs.open()))
    series = list(csv.DictReader(SERIES.open()))
    assert len(runs) == len(RECS) * 28 * len(SEEDS) * 2
    scores = {
        (r["receptor"], r["compound"], int(r["seed"]), r["arm"]): (
            float(r["top_score"]) if r["status"] == "0" and r.get("top_score") else None
        )
        for r in runs
    }
    assert len(scores) == len(runs)
    compounds = []
    for rec, r in itertools.product(RECS, series):
        y = 9 - math.log10(float(r["ki_nm"]))
        row = dict(
            receptor=rec,
            compound=r["ligand"],
            activity_id=r["activity"],
            ki_nm=r["ki_nm"],
            pki=y,
            heavy=int(r["heavy"]),
            charge=int(r["charge"]),
            prepared=int(r["included"]),
        )
        for arm in ("baseline", "fixed"):
            cells = [scores.get((rec, r["ligand"], s, arm)) for s in SEEDS]
            row[arm + "_completed_seeds"] = sum(v is not None for v in cells)
            # Exact arithmetic to half the printed score resolution; avoid floating-point false ties.
            row[arm + "_strength"] = (
                round(-float(np.mean(cells)), 4)
                if all(v is not None for v in cells)
                else None
            )
        compounds.append(row)
    write(args.out / "compounds.csv", compounds)
    pairs, metrics, point, draws = [], [], [], []
    for rec in RECS:
        local = [r for r in compounds if r["receptor"] == rec]
        finite = [
            r
            for r in local
            if all(r[a + "_strength"] is not None for a in ("baseline", "fixed"))
        ]
        local_pairs = []
        for u, v in itertools.combinations(finite, 2):
            dy = u["pki"] - v["pki"]
            strong, weak = (u, v) if dy >= 0 else (v, u)
            fold = float(weak["ki_nm"]) / float(strong["ki_nm"])
            pr = dict(
                receptor=rec,
                compound1=strong["compound"],
                compound2=weak["compound"],
                fold=fold,
                same_charge_size=int(
                    strong["heavy"] == weak["heavy"]
                    and strong["charge"] == weak["charge"]
                ),
            )
            for arm in ("baseline", "fixed"):
                delta = round(strong[arm + "_strength"] - weak[arm + "_strength"], 4)
                pr[arm + "_gap"] = delta
                pr[arm + "_credit"] = (
                    None if dy == 0 else (1 if delta > 0 else 0.5 if delta == 0 else 0)
                )
            local_pairs.append(pr)
        pairs.extend(local_pairs)
        y = np.asarray([r["pki"] for r in finite])
        x = np.asarray(
            [[r[a + "_strength"] for r in finite] for a in ("baseline", "fixed")]
        )
        est = np.asarray([ranking(y, v) for v in x])
        for arm, values in zip(("baseline", "fixed"), est):
            row = dict(
                receptor=rec,
                arm=arm,
                source_n=len(local),
                expected_prepared_n=28,
                paired_n=len(finite),
                complete_n=sum(r[arm + "_strength"] is not None for r in local),
                rho=values[0],
                tau_b=values[1],
                all_pairs=len(local_pairs),
                experimental_ties=sum(r["fold"] == 1 for r in local_pairs),
            )
            for subset, threshold, matched in (
                ("non_ties", 1, False),
                ("fold3", 3, False),
                ("matched_fold3", 3, True),
            ):
                eligible = [
                    r
                    for r in local_pairs
                    if r["fold"] > 1
                    and r["fold"] >= threshold
                    and (not matched or r["same_charge_size"])
                ]
                row[subset + "_n"] = len(eligible)
                row[subset + "_credit"] = sum(r[arm + "_credit"] for r in eligible)
                row[subset + "_prediction_ties"] = sum(
                    r[arm + "_gap"] == 0 for r in eligible
                )
            metrics.append(row)
        # Both receptor rows have identical compounds. The identical RNG stream resamples
        # complete compound blocks jointly across both receptors, seeds and arms.
        rng = np.random.default_rng(173092)
        samples = np.array(
            [rng.integers(0, len(y), len(y)) for _ in range(args.resamples)]
        )
        boot = np.array(
            [ranking(y[idx], x[1, idx]) - ranking(y[idx], x[0, idx]) for idx in samples]
        )
        point.append(est[1] - est[0])
        draws.append(boot)
    write(args.out / "pairs.csv", pairs)
    write(args.out / "metrics.csv", metrics)
    intervals = []
    for group, values, effects in zip(
        (*RECS, "equal_receptor_macro"),
        (*draws, np.mean(draws, axis=0)),
        (*point, np.mean(point, axis=0)),
    ):
        for j, label in enumerate(("rho", "tau_b", "fold3_concordance")):
            d = values[:, j]
            d = d[np.isfinite(d)]
            lo, hi = np.quantile(d, [0.025, 0.975])
            intervals.append(
                dict(
                    group=group,
                    metric=label,
                    fixed_minus_baseline=effects[j],
                    ci025=lo,
                    ci975=hi,
                    valid_bootstrap=len(d),
                    requested_bootstrap=args.resamples,
                )
            )
    write(args.out / "bootstrap.csv", intervals)


if __name__ == "__main__":
    main()
