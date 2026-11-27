"""Within-series uncertainty for paired per-compound changes in >=3-fold Kd concordance.

Draw N ligand identities with replacement from the frozen N=47, keeping both
receptors and all score arms together. An eligible *distinct* i,j comparison
has multiplicity n_i*n_j in that draw. Tie credit is one half; each replicate
divides by its own eligible multiplicity. The percentile interval describes
sampling compounds from this one series, conditional on the assay, structures,
single seed and fixed scoring protocol. It excludes assay measurement error,
scaffold-, target-, seed- and receptor-population uncertainty. It is not a
pairwise-binomial interval. A replicate without eligible pairs is omitted and
counted. Receptors and consensus share draws; their rows are dependent.

Run: python compound_bootstrap.py [directory containing manifest.csv, scores.csv
and summary.csv] [--out destination.csv]. No new docking or extra data needed.
"""

import argparse
import csv
import itertools
from decimal import Decimal
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RECEIVERS = ("1KJL", "3T1L")
SEED = 20261127
DRAWS = 10000


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def first_scores(folder):
    values = {}
    methods = set()
    for row in read(folder / "scores.csv"):
        if row["receiver"] not in RECEIVERS or row["method"] == "size":
            continue
        key = (row["method"], row["receiver"], int(row["compound"]))
        ticks = Decimal(row["score"]) * 1000
        if key in values or ticks != int(ticks):
            raise ValueError(f"Duplicate or non-PDBQT precision score for {key}")
        methods.add(key[0])
        values[key] = int(ticks)
    return values, sorted(methods - {"parent"})


def run(folder):
    manifest = sorted(read(folder / "manifest.csv"), key=lambda r: int(r["compound"]))
    ids = [int(r["compound"]) for r in manifest]
    n = len(ids)
    if len(set(ids)) != n or len({r["smiles"] for r in manifest}) != n or n != 47:
        raise ValueError("Expected the frozen panel of 47 distinct graphs")
    kd = [Decimal(r["kd_uM"]) for r in manifest]
    size_charge = [(int(r["heavy"]), int(r["charge"])) for r in manifest]
    scores, methods = first_scores(folder)
    if len(scores) != n * 2 * (len(methods) + 1):
        raise ValueError("Incomplete paired scores")
    observed = {
        (r["receiver"], r["method"], r["subset"]): (
            int(r["n_pairs"]),
            float(r["credit"]),
        )
        for r in read(folder / "summary.csv")
    }
    # Draw same vector of multiplicities for all receiver, arm and subset comparisons.
    rng = np.random.default_rng(SEED)
    counts = rng.multinomial(n, np.repeat(1 / n, n), size=DRAWS)
    results = []
    for subset in ("all", "equal_heavy_charge"):
        edges = [
            (i, j)
            for i, j in itertools.permutations(range(n), 2)
            if kd[j] >= 3 * kd[i]
            and (subset == "all" or size_charge[i] == size_charge[j])
        ]
        strong, weak = np.array(edges).T
        weights = counts[:, strong] * counts[:, weak]
        den = weights.sum(axis=1)
        valid = den > 0
        for receiver in (*RECEIVERS, "consensus"):
            recs = RECEIVERS if receiver == "consensus" else (receiver,)
            credits = {}
            for method in ["parent", *methods]:
                # The sum preserves exact ranking under the 2-receiver arithmetic mean.
                score = np.array(
                    [sum(scores[method, r, idx] for r in recs) for idx in ids]
                )
                credits[method] = 1 + np.sign(score[weak] - score[strong])
                if observed[receiver, method, subset] != (
                    len(edges),
                    credits[method].sum() / 2,
                ):
                    raise ValueError(
                        f"Pair reconstruction differs for {receiver}/{method}/{subset}"
                    )
            for arm in methods:
                delta = credits[arm] - credits["parent"]
                # A signed half-credit difference over multiplicities, in percentage points.
                changes = 50 * (weights @ delta)[valid] / den[valid]
                lo, med, hi = np.percentile(changes, (2.5, 50, 97.5))
                results.append(
                    dict(
                        receiver=receiver,
                        arm_minus_parent=arm,
                        subset=subset,
                        panel_ligands=n,
                        ligands_on_eligible_edges=len(set(strong) | set(weak)),
                        eligible_edges=len(edges),
                        parent_credit=credits["parent"].sum() / 2,
                        arm_credit=credits[arm].sum() / 2,
                        observed_delta_pp=round(50 * delta.mean(), 4),
                        lower_2p5_pp=round(lo, 4),
                        median_pp=round(med, 4),
                        upper_97p5_pp=round(hi, 4),
                        resampled_ligands_per_draw=n,
                        draws=DRAWS,
                        valid_draws=int(valid.sum()),
                        seed=SEED,
                    )
                )
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", nargs="?", type=Path, default=HERE)
    parser.add_argument("--out", type=Path, default=HERE / "compound_bootstrap.csv")
    args = parser.parse_args()
    rows = run(args.folder)
    with args.out.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
