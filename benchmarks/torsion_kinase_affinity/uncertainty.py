"""Within-target rank correlations and paired whole-compound bootstrap."""

import argparse
import csv
import itertools
from pathlib import Path
from statistics import median

import numpy as np
from scipy.stats import kendalltau, spearmanr

from run import ARMS, STATES


def main(work):
    panel = list(csv.DictReader((work / "panel.csv").open()))
    cells = list(csv.DictReader((work / "cells.csv").open()))
    out = []
    for target, states in STATES.items():
        rows = [r for r in panel if r["target"] == target]
        kd = np.array([float(r["kd_nm"]) for r in rows])
        pairs = [
            (i, j)
            for i, j in itertools.permutations(range(len(rows)), 2)
            if kd[j] >= 3 * kd[i]
        ]
        i, j = np.array(pairs).T
        genotype = np.random.default_rng(61309).multinomial(
            len(rows), np.ones(len(rows)) / len(rows), 5000
        )
        weight = genotype[:, i] * genotype[:, j]
        weight = weight[weight.sum(axis=1) > 0]
        for arm in ARMS:
            e = np.array(
                [
                    np.mean(
                        [
                            median(
                                float(c["pool_score"])
                                for c in cells
                                if c["cid"] == r["cid"]
                                and c["state"] == s
                                and c["arm"] == arm
                            )
                            for s in states
                        ]
                    )
                    for r in rows
                ]
            )
            diff = e[j] - e[i]
            concordance = np.where(np.abs(diff) < 1e-8, 0.5, diff > 0)
            if arm == "native":
                reference = concordance
            estimates = (weight * (concordance - reference)).sum(axis=1) / weight.sum(
                axis=1
            )
            lo, hi = np.quantile(estimates, [0.025, 0.975])
            out.append(
                dict(
                    target=target,
                    arm=arm,
                    spearman=spearmanr(kd, e).statistic,
                    kendall_tau_b=kendalltau(kd, e).statistic,
                    threefold_accuracy=concordance.mean(),
                    paired_delta_pp=100 * (concordance - reference).mean(),
                    bootstrap_low_pp=100 * lo,
                    bootstrap_high_pp=100 * hi,
                    draws=len(estimates),
                )
            )
    with (work / "uncertainty.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=out[0])
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    main(p.parse_args().work)
