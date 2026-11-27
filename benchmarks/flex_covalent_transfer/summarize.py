"""Recover continuous distributions, per-pair effects, and within-assay ranks from archived rows."""

import argparse
import csv
import itertools
import math
from pathlib import Path
import statistics

import numpy as np
from scipy.stats import kendalltau, spearmanr, rankdata

from experiment import csvout, SEEDS
from prepare import MEASUREMENTS


def dist(values):
    if not values:
        return dict(n=0)
    return dict(
        n=len(values),
        mean=float(np.mean(values)),
        **dict(
            zip(
                ("min", "q25", "median", "q75", "p90", "p95", "max"),
                np.quantile(values, [0, 0.25, 0.5, 0.75, 0.9, 0.95, 1]).tolist(),
            )
        ),
    )


def structure(out, records):
    cases = [r for r in records if r["cohort"] != "assay"]

    def groups(row):
        return [row["cohort"] + "/all", row["cohort"] + "/" + row["target"]]

    stats, ecdf, effects, raw = [], [], [], []
    for grp in sorted({g for r in cases for g in groups(r)}):
        rows = [r for r in cases if grp in groups(r)]
        for arm in ("baseline", "fixed", "official"):
            sub = [r for r in rows if r["arm"] == arm]
            for k in (
                "top1",
                "best5",
                "best_reported",
                "regret",
                "reported",
                "wall_s",
                "cpu_s",
                "peak_kib",
                "hunt_evals",
                "second_evals",
                "explicit_evals",
            ):
                values = [float(r[k]) for r in sub if r.get(k)]
                stats.append(
                    dict(
                        group=grp,
                        arm=arm,
                        quantity=k,
                        **dist(values),
                        denominator=len(sub),
                    )
                )
                if k in ("top1", "best5", "best_reported", "regret"):
                    for i, v in enumerate(sorted(values), 1):
                        ecdf.append(
                            dict(
                                group=grp,
                                arm=arm,
                                quantity=k,
                                value=v,
                                f=i / len(values),
                            )
                        )
        pairkey = lambda r: (r["query"], r["receiver"], int(r["seed"]))
        arms = {
            a: {pairkey(r): r for r in rows if r["arm"] == a}
            for a in ("baseline", "fixed")
        }
        for k in ("top1", "best5", "best_reported", "regret"):
            diffs = []
            for ident, a in arms["baseline"].items():
                b = arms["fixed"].get(ident, {})
                if not a.get(k) or not b.get(k):
                    continue
                d = float(b[k]) - float(a[k])
                diffs.append((ident[2], d))
                raw.append(
                    dict(
                        group=grp,
                        query=ident[0],
                        receiver=ident[1],
                        seed=ident[2],
                        quantity=k,
                        baseline=a[k],
                        fixed=b[k],
                        delta=d,
                    )
                )
            byseed = {s: [d for seed, d in diffs if seed == s] for s in SEEDS}
            boot = [
                np.mean([v for seed in seeds for v in byseed[seed]])
                for seeds in itertools.product(SEEDS, repeat=len(SEEDS))
            ]
            lo, hi = np.quantile(boot, [0.025, 0.975])
            effects.append(
                dict(
                    group=grp,
                    quantity=k,
                    **dist([d for s, d in diffs]),
                    smaller=sum(d < 0 for s, d in diffs),
                    larger=sum(d > 0 for s, d in diffs),
                    zero=sum(d == 0 for s, d in diffs),
                    seed_block_mean_low=lo,
                    seed_block_mean_high=hi,
                )
            )
    for n, data in (
        ("distribution", stats),
        ("ecdf", ecdf),
        ("effects", effects),
        ("paired", raw),
    ):
        csvout(out / (n + ".csv"), data)


def affinity(out, records):
    measured = list(csv.DictReader(MEASUREMENTS.open()))
    data = [r for r in records if r["cohort"] == "assay"]
    metrics, pairs, scores = [], [], []
    for rx, arm in itertools.product(
        ("2NNQ", "3FR2", "equal_mean"), ("baseline", "fixed")
    ):
        vals = {}
        for m in measured:
            rows = [
                r
                for r in data
                if r["arm"] == arm
                and r["query"] == m["id"]
                and (r["receiver"] == rx or rx == "equal_mean")
                and r["top1_score"]
            ]
            if len(rows) == (6 if rx == "equal_mean" else 3):
                vals[m["id"]] = statistics.mean(float(r["top1_score"]) for r in rows)
            scores.append(
                dict(
                    receiver=rx,
                    arm=arm,
                    id=m["id"],
                    kd_uM=m["Kd_uM"],
                    pKd=6 - math.log10(float(m["Kd_uM"])),
                    mean_score=vals.get(m["id"]),
                    n=len(rows),
                )
            )
        eligible = []
        for a, b in itertools.combinations(measured, 2):
            ka, kb = float(a["Kd_uM"]), float(b["Kd_uM"])
            use = abs(ka - kb) > float(a["reported_error_uM"]) + float(
                b["reported_error_uM"]
            )
            da = (
                vals[b["id"]] - vals[a["id"]]
                if {a["id"], b["id"]} <= vals.keys()
                else None
            )
            c = (
                (float(da * (kb - ka) > 0) if da else 0.5)
                if use and da is not None
                else None
            )
            if c is not None:
                eligible.append(c)
            pairs.append(
                dict(
                    receiver=rx,
                    arm=arm,
                    a=a["id"],
                    b=b["id"],
                    fold=max(ka / kb, kb / ka),
                    eligible=use,
                    score_delta=da,
                    agreement=c,
                )
            )
        aa = [
            r
            for r in scores
            if r["receiver"] == rx and r["arm"] == arm and r["mean_score"] is not None
        ]
        obs = [r["pKd"] for r in aa]
        pred = [-r["mean_score"] for r in aa]
        metrics.append(
            dict(
                receiver=rx,
                arm=arm,
                n=len(aa),
                planned=len(measured),
                spearman=spearmanr(obs, pred).statistic,
                kendall_b=kendalltau(obs, pred).statistic,
                concordant=sum(eligible),
                pairs=len(eligible),
                all_pairs=math.comb(len(measured), 2),
            )
        )
    intervals = []
    for rx in ("2NNQ", "3FR2", "equal_mean"):
        arrays = {
            a: [r for r in scores if r["receiver"] == rx and r["arm"] == a]
            for a in ("baseline", "fixed")
        }
        observed = np.array([r["pKd"] for r in arrays["baseline"]])
        predictions = {
            a: -np.array([r["mean_score"] for r in rr]) for a, rr in arrays.items()
        }
        # Exhaust all ligand-cluster bootstrap samples. Both methods resample identical compounds.
        samples = np.array(
            list(itertools.product(range(len(observed)), repeat=len(observed)))
        )
        samples = samples[np.ptp(observed[samples], axis=1) > 0]

        def normalized_rank(v):
            rv = rankdata(v[samples], axis=1)
            centered = rv - rv.mean(axis=1, keepdims=True)
            return centered / np.linalg.norm(centered, axis=1, keepdims=True)

        ranked = normalized_rank(observed)
        delta = np.sum(
            ranked
            * (
                normalized_rank(predictions["fixed"])
                - normalized_rank(predictions["baseline"])
            ),
            axis=1,
        )
        lo, hi = np.quantile(delta, [0.025, 0.975])
        intervals.append(
            dict(
                receiver=rx,
                quantity="delta_spearman",
                draws=len(delta),
                all_draws=len(observed) ** len(observed),
                percentile_025=lo,
                percentile_975=hi,
            )
        )
    csvout(out / "assay_intervals.csv", intervals)
    for n, data in (
        ("assay_metrics", metrics),
        ("assay_pairs", pairs),
        ("assay_compounds", scores),
    ):
        csvout(out / (n + ".csv"), data)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    w = p.parse_args().out
    rows = list(csv.DictReader((w / "cases.csv").open()))
    structure(w, rows)
    if any(r["cohort"] == "assay" for r in rows):
        affinity(w, rows)


if __name__ == "__main__":
    main()
