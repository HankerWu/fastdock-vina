"""Paired molecule bootstrap and leave-one-molecule-out influence on Gal3."""

import argparse
import csv
from pathlib import Path
import numpy as np

BENCH = Path(__file__).resolve().parents[1]


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def audit(work):
    manifest = read(BENCH / "pair_gradient/galectin_external/manifest.csv")
    ids = [r["compound"] for r in manifest]
    idx = {x: i for i, x in enumerate(ids)}
    groups = {}
    for r in read(work / "pairs.csv"):
        if r["receptor"] == "mean" and r["reliable"] == "True":
            groups.setdefault(r["arm"], {})[idx[r["strong"]], idx[r["weak"]]] = float(
                r["credit"]
            )
    edges = list(groups["parent"])
    i, j = np.array(edges).T
    size = np.array([int(r["heavy"]) for r in manifest])
    size_gap = size[i] - size[j]
    groups["size"] = dict(
        zip(edges, (size_gap > 0).astype(float) + 0.5 * (size_gap == 0))
    )
    weights = np.random.default_rng(623113).multinomial(
        len(ids), [1 / len(ids)] * len(ids), size=10000
    )
    pair_weights = weights[:, i] * weights[:, j]
    parent = np.array([groups["parent"][edge] for edge in edges])
    out = []
    for arm, group in groups.items():
        scores = np.array([group[edge] for edge in edges])
        diff = scores - parent
        fractions = (pair_weights @ diff) / pair_weights.sum(axis=1)
        lo, hi = np.quantile(fractions, [0.025, 0.975])
        jackknife = [diff[(i != k) & (j != k)].sum() for k in range(len(ids))]
        keep = np.array([key not in {"45", "46"} for key in ids])
        subset = keep[i] & keep[j]
        out.append(
            dict(
                arm=arm,
                pairs=len(edges),
                credit=scores.sum(),
                delta=diff.sum(),
                paired_fraction_low=lo,
                paired_fraction_high=hi,
                min_delete_one=min(jackknife),
                max_delete_one=max(jackknife),
                without_bivalent_pairs=subset.sum(),
                without_bivalent_credit=scores[subset].sum(),
            )
        )
    with (work / "influence.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(out[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    audit(p.parse_args().work)
