"""Assay-level ranks, paired factorial contrasts and ligand-clustered confidence intervals."""
import argparse
import itertools
from collections import defaultdict
from decimal import Decimal
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau, spearmanr
from replay import load, save

ARMS = ("native", "strain", "rotor", "joint", "strain_unscaled_no_rotor")


def group(rows, fields):
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[k] for k in fields)].append(row)
    return groups


def transfer_coverage(raw):
    coverage = []
    for (target,), rows in sorted(group(raw, ["target"]).items()):
        compounds = group(rows, ["compound"])
        values = {k: Decimal(v[0]["measurement"]) for k, v in compounds.items()}
        complete = {k for k, v in compounds.items() if all(r["rotor_status"] == "ok" for r in v)}
        edges = [(a, b) for a, b in itertools.combinations(values, 2)
                 if max(values[a], values[b]) >= 3 * min(values[a], values[b])]
        coverage.append(dict(target=target, compounds=len(compounds), common_compounds=len(complete),
                             poses=len(rows), common_poses=sum(len(compounds[k]) for k in complete),
                             reliable_pairs=len(edges), common_reliable_pairs=sum(a in complete and b in complete for a,b in edges),
                             abstaining_compounds=" ".join(sorted(k[0] for k in compounds if k not in complete))))
    return coverage


def analyze(rows, name, estimator, boot):
    by_ligand = group(rows, ["compound"])
    labels = sorted(k[0] for k in by_ligand)
    measurements, descriptors, chem = [], [], []
    blocks = []
    for label in labels:
        ordered = sorted(by_ligand[label,], key=lambda r: (r["receptor"], r["seed"]))
        fixed = {(Decimal(str(r["measured_nM"])), r["heavy"], r["charge"]) for r in ordered}
        if len(fixed) != 1:
            raise ValueError("Conflicting identities or labels")
        nm, heavy, charge = fixed.pop()
        measurements.append(nm)
        chem.append((int(heavy), int(charge)))
        blocks.append([(r["receptor"], r["seed"]) for r in ordered])
        descriptors.append([[float(r[a]) for a in ARMS] for r in ordered])
    if any(b != blocks[0] for b in blocks):
        raise ValueError("Ligands have different receptor/seed blocks")
    if len(set(blocks[0])) != len(blocks[0]):
        raise ValueError("Duplicate ligand cases")
    raw = np.array(descriptors)
    receptors = np.array([r for r, s in blocks[0]])
    scores = np.mean([np.median(raw[:, receptors == rec], axis=1) for rec in sorted(set(receptors))], axis=0)
    high, low = [], []
    for i, j in itertools.combinations(range(len(labels)), 2):
        if measurements[i] == measurements[j]:
            continue
        strong, weak = sorted((i, j), key=lambda k: measurements[k])
        high.append(strong)
        low.append(weak)
    strong, weak = np.array(high), np.array(low)
    gaps = (np.median(raw[weak] - raw[strong], axis=1) if estimator == "paired_gap_median"
            else scores[weak] - scores[strong])
    credits = np.where(abs(gaps) < 1e-10, .5, (gaps > 0).astype(float))
    edges = []
    for pos, (i, j) in enumerate(zip(strong, weak)):
        reliable = measurements[j] >= 3 * measurements[i]
        same = chem[i] == chem[j]
        edges.append(dict(**name, estimator=estimator, strong=labels[i], weak=labels[j],
                          fold=float(measurements[j]/measurements[i]), threefold=int(reliable),
                          same_heavy_charge=int(same), **dict(zip(ARMS, gaps[pos]))))
    keep = np.array([bool(e["threefold"]) for e in edges])
    matched = np.array([bool(e["same_heavy_charge"]) for e in edges]) & keep
    if not keep.any():
        raise ValueError("Panel has no threefold comparisons")
    counts = credits[keep].sum(axis=0)
    matched_credit = credits[matched].sum(axis=0)
    result = dict(**name, estimator=estimator, n_compounds=len(labels), n_nontied=len(edges),
                  n_3fold=int(sum(keep)), n_3fold_size_charge=int(sum(matched)))
    for k, arm in enumerate(ARMS):
        result[arm + "_credit"] = counts[k]
        result[arm + "_size_charge_credit"] = matched_credit[k]
        scalar = estimator != "paired_gap_median"
        x = np.array([float(v) for v in measurements])
        result[arm + "_rho"] = float(spearmanr(x, scores[:,k]).statistic) if scalar else ""
        result[arm + "_tau_b"] = float(kendalltau(x, scores[:,k]).statistic) if scalar else ""
    interaction = credits[:,3] - credits[:,1] - credits[:,2] + credits[:,0]
    contrasts = np.c_[credits[:,3] - credits[:,1], credits[:,3] - credits[:,2], interaction]
    names = ("joint_minus_strain", "joint_minus_rotor", "interaction")
    # Draw ligand identities; a pair's multiplicity is the product of its endpoints' multiplicities.
    # The same draw is applied to all arms and their paired receptor and seed observations.
    draws = np.random.default_rng(241127).multinomial(len(labels), [1/len(labels)] * len(labels), size=boot)
    weights = draws[:, strong[keep]] * draws[:, weak[keep]]
    totals = weights.sum(axis=1)
    valid = totals > 0
    samples = weights[valid] @ contrasts[keep] / totals[valid, None]
    result["bootstrap_n_valid"] = int(sum(valid))
    for k, contrast in enumerate(names):
        result[contrast + "_credit"] = contrasts[keep,k].sum()
        result[contrast + "_pp"] = 100 * contrasts[keep,k].mean()
        lo, hi = 100 * np.quantile(samples[:,k], [.025, .975])
        result[contrast + "_ci_low_pp"] = lo
        result[contrast + "_ci_high_pp"] = hi
    return result, edges


def main(out, boot):
    save(out / "transfer_coverage.csv", transfer_coverage(load(out / "transfer.csv")))
    summary, pairs = [], []
    for key, rows in sorted(group(load(out / "first_poses.csv.xz"), ["panel", "target", "assay", "endpoint"]).items()):
        common = dict(zip(("panel", "target", "assay", "endpoint"), key))
        by_rec = group(rows, ["receptor"])
        groups = [("consensus", rows, "mean_receptor_score")]
        groups.extend((key[0], sub, "single_receptor") for key, sub in sorted(by_rec.items()))
        if common["target"] == "FXa":
            groups.append(("both", rows, "paired_gap_median"))
        for rec, sub, estimator in groups:
            result, edges = analyze(sub, dict(**common, receptor=rec), estimator, boot)
            summary.append(result)
            pairs.extend(edges)
    save(out / "summary.csv", summary)
    save(out / "pair_gaps.csv.xz", pairs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--bootstrap", type=int, default=5000)
    args = p.parse_args()
    main(args.out, args.bootstrap)
