"""Can one two-pose oracle per receptor realize all individually favorable choices?"""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog

from replay import KEY, PRESENT, RT, TOL, indicators, partition_kernels, recover_factors


def main(source, pairs_path, output):
    sites = pd.read_csv(source, float_precision="round_trip")
    pairs = pd.read_csv(pairs_path)
    rows, weights = [], []
    for receiver, pairs_rec in pairs.query("result == 'reversible'").groupby("receiver"):
        # All poses see the same dry receptor and the same labelled set of sites.
        group = sites[sites.receiver == receiver]
        names = sorted(group.site.unique())
        pairs_a, gaps, preferred = [], [], []
        ids = []
        for row in pairs_rec.itertuples(index=False):
            part = group
            for col in KEY:
                part = part[part[col] == getattr(row, col)]
            u = part.set_index("site").loc[names, ["penalty_original","penalty_winner"]].to_numpy().T
            pairs_a.append(recover_factors(u))
            gaps.append(row.dry_gap)
            preferred.append(1 if row.rmsd_change < 0 else -1)
            ids.append({k: getattr(row, k) for k in KEY})
        aa = np.asarray(pairs_a)
        loss = np.max(1-aa, axis=(0,1))
        ids_nonzero = np.flatnonzero(loss > 0)
        by_impact = ids_nonzero[np.argsort(-loss[ids_nonzero], kind="stable")]
        used, others = by_impact[:12], by_impact[12:]
        bits = indicators(len(used))
        # State columns model all other sites independently, providing an
        # explicit full-site distribution (a subset of all admissible laws).
        contrast = []
        for a, gap, sign in zip(aa, gaps, preferred, strict=True):
            k = [partition_kernels(bits, row[used]) *
                 np.prod(1-PRESENT*(1-row[others])) for row in a]
            contrast.append(sign * (np.exp(-gap/RT)*k[1] - k[0]))
        contrast = np.stack(contrast)
        eq = np.vstack((np.ones(len(bits)), bits.T))
        target = np.r_[1.,np.full(len(used), PRESENT)]
        # Each missing site's marginal bounds the possible change in expected
        # insertion weight between two extensions of the same retained law.
        residual = PRESENT * (np.sum(1-aa[:, 0, others], axis=1)
                              + np.exp(-np.array(gaps)/RT) * np.sum(1-aa[:, 1, others], axis=1))
        eq_with_t = np.column_stack((eq,np.zeros(len(eq))))
        best_count = -1
        upper_count = -1
        q_best = None
        exclusion_bound = None
        for count in range(len(aa), -1, -1):
            for subset in itertools.combinations(range(len(aa)), count):
                if not subset:
                    raise ArithmeticError("Even a singleton produced no positive witness")
                chosen = list(subset)
                sub = contrast[chosen]
                r = linprog(np.r_[np.zeros(len(bits)), -1.0],
                            A_ub=np.column_stack((-sub, np.ones(count))),
                            b_ub=np.zeros(count), A_eq=eq_with_t, b_eq=target,
                            bounds=[(0,None)]*len(bits)+[(None,None)], method="highs")
                if not r.success:
                    raise RuntimeError(r.message)
                q = r.x[:-1]
                if np.max(np.abs(eq @ q-target)) > 1e-9 or np.min(q) < -1e-10:
                    raise ArithmeticError("Infeasible common occupation witness")
                signed = contrast @ q
                if (signed[chosen] > TOL).all():
                    best_count = int((signed > TOL).sum())
                    upper_count = max(upper_count, count)
                    q_best = q
                    break
                # Maximum minimum margin cannot exceed its value averaged with
                # any probability vector across the selected constraints.
                w = np.maximum(-r.ineqlin.marginals, 0)
                w /= w.sum()
                h = w @ sub
                y = -r.eqlin.marginals
                upper = (target @ y + max(0.0, np.max(h - eq.T @ y))
                         + w @ residual[chosen])
                if upper >= -TOL:
                    upper_count = max(upper_count, count)
                if count == len(aa):
                    exclusion_bound = float(upper)
            if q_best is not None:
                break
        q = q_best
        margins = contrast @ q
        for key, margin in zip(ids, margins, strict=True):
            rows.append(dict(**key, desired_margin=margin))
        weights.append(dict(receiver=receiver, correlated_sites=[names[i] for i in used],
                            other_sites="independent p=0.75",
                            support=[dict(occupied=[names[used[j]] for j in np.flatnonzero(bits[i])],
                                          probability=q[i]) for i in np.flatnonzero(q>0)],
                            reversible_pairs=len(aa), largest_demonstrated_count=best_count,
                            upper_count=upper_count, full_set_upper_margin=exclusion_bound))
    pd.DataFrame(rows).to_csv(output / "shared_margins.csv", index=False, float_format="%.12g")
    (output / "shared_witnesses.json").write_text(json.dumps(weights, indent=2)+"\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("pairs", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    main(args.source, args.pairs, args.output)
