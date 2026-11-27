"""Bound how a common distribution of occupied water sites can rank two frozen poses."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog

# The same geometric site model and temperature as benchmarks/mobile_water.
PRESENT = 0.75
RT = 0.00198720425864083 * 298.15
KEY = ["panel", "target", "cap", "query", "receiver", "seed", "search"]
TOL = 1e-8


def indicators(n):
    codes = np.arange(2**n, dtype=np.uint32)[:, None]
    return ((codes >> np.arange(n, dtype=np.uint32)) & 1).astype(float)


def partition_kernels(states, availability):
    kernels = np.ones(len(states))
    for occupied, available in zip(states.T, availability, strict=True):
        kernels *= 1.0 - occupied * (1.0 - available)
    return kernels


def lower_bound(availability, gap, keep=12):
    """Positive min margin certifies the water-selected pose for every joint law.

    The primal gives a feasible independent extension of the retained sites.
    The dual is checked against every retained occupation pattern. The one-sided
    union bound on discarded sites remains valid for arbitrary dependence.
    """
    if not np.all((0 <= availability) & (availability <= 1)):
        raise ValueError("Invalid insertion factors")
    scale = np.exp(-gap / RT)
    deficit = 1.0 - availability
    order = np.argsort(-np.max(deficit, axis=0), kind="stable")
    nontrivial = order[np.max(deficit[:, order], axis=0) > 0.0]
    used, dropped = nontrivial[:keep], nontrivial[keep:]
    bits = indicators(len(used))
    A = np.vstack((np.ones(len(bits)), bits.T))
    b = np.r_[1.0, np.full(len(used), PRESENT)]
    ko, kw = [partition_kernels(bits, row[used]) for row in availability]
    c = scale * kw - ko
    result = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method="highs",
                     options={"primal_feasibility_tolerance": 1e-10,
                              "dual_feasibility_tolerance": 1e-10})
    if not result.success:
        raise RuntimeError(result.message)
    q = result.x
    dual = result.eqlin.marginals
    if np.max(np.abs(A @ q - b)) > 1e-9 or np.min(q) < -1e-10:
        raise ArithmeticError("Infeasible joint occupation witness")
    dual_violation = max(0.0, np.max(A.T @ dual - c))
    # q is shared by both poses. Bind omitted sites independently only for
    # this witness; the lower bound covers all couplings to the omitted sites.
    extension = np.prod(1.0 - PRESENT * deficit[:, dropped], axis=1)
    witness = scale * (kw @ q) * extension[1] - (ko @ q) * extension[0]
    truncation = scale * PRESENT * deficit[1, dropped].sum()
    lower = b @ dual - dual_violation - truncation
    if lower > witness + 1e-9:
        raise ArithmeticError("Inverted primal-dual bounds")
    return dict(active_states=len(c), projected_sites=len(used), lower=lower,
                feasible_margin=witness, lower_truncation=truncation,
                dual_violation=dual_violation)


def recover_factors(penalties):
    # exp(-penalty/RT) = 1-p+p*a. expm1 preserves weak interactions.
    values = 1.0 + np.expm1(-penalties / RT) / PRESENT
    if np.min(values) < -1e-12 or np.max(values) > 1.0 + 1e-12:
        raise ValueError("Historical penalty outside its single-site range")
    return np.clip(values, 0, 1)


def numeric_checks():
    # Two sites with equal marginal p have one free parameter t=P(11).
    # For m=2 any linear objective attains its minimum at either endpoint.
    p = PRESENT
    for avail in (np.array([[.2, .3], [.8, .4]]), np.array([[0., .9], [.2, 0.]])):
        scale = np.exp(-.13 / RT)
        kernels = [partition_kernels(indicators(2), row) for row in avail]
        c = scale * kernels[1] - kernels[0]
        expected = min(np.array([1 - 2*p + t, p-t, p-t, t]) @ c
                       for t in [max(0, 2*p-1), p])
        r = lower_bound(avail, .13)
        if max(abs(r[k] - expected) for k in ["lower", "feasible_margin"]) > 1e-9:
            raise ArithmeticError("Two-site Fréchet bounds disagree")
        approx = lower_bound(avail, .13, keep=1)
        if not approx["lower"] - TOL <= expected <= approx["feasible_margin"] + TOL:
            raise ArithmeticError("Truncation excludes the exact optimum")


def main(source, out):
    numeric_checks()
    sites_file = source / "changed_winner_sites.csv"
    cells_file = source / "cells.csv"
    groups = pd.read_csv(sites_file, float_precision="round_trip").groupby(KEY)
    cells = pd.read_csv(cells_file).set_index(KEY + ["ranker"])
    rows = []
    for key, group in groups:
        gap, = group.dry_gap.unique()
        winner, = group.winner.unique()
        original = cells.loc[key + ("dry",)]
        selected = cells.loc[key + ("water",)]
        if original.model != 1 or selected.model != winner or len(group) != original.n_sites:
            raise ValueError(f"Wrong pair or lost sites: {key}")
        u = group[["penalty_original", "penalty_winner"]].to_numpy().T
        a = recover_factors(u)
        ind = np.exp(-gap / RT) * np.exp(-u[1].sum()/RT) - np.exp(-u[0].sum()/RT)
        if ind <= 0:
            raise ArithmeticError(f"Original switch not reproduced: {key}")
        result = lower_bound(a, gap)
        conclusion = ("forced" if result["lower"] > TOL else
                      "reversible" if result["feasible_margin"] < -TOL else "unresolved")
        rows.append(dict(zip(KEY, key), original_model=1, water_model=winner,
                         n_sites=len(group), dry_gap=gap, original_rmsd=original.rmsd,
                         water_rmsd=selected.rmsd, rmsd_change=selected.rmsd-original.rmsd,
                         independent_margin=ind, **result, result=conclusion))
    frame = pd.DataFrame(rows)
    expected = cells.reset_index().query("ranker == 'water' and model != 1")
    if set(groups.groups) != set(expected.set_index(KEY).index):
        raise ValueError("The compared pairs do not cover all changed winners")
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out / "pairs.csv", index=False, float_format="%.12g")
    info = {"input_sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                              for f in (sites_file, cells_file)},
            "pairs": len(frame), "forced": int((frame.result == "forced").sum()),
            "reversible": int((frame.result == "reversible").sum()),
            "unresolved": int((frame.result == "unresolved").sum())}
    (out / "summary.json").write_text(json.dumps(info, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path, help="benchmarks/soft_water_interaction/results")
    p.add_argument("out", type=Path)
    args = p.parse_args()
    main(args.source, args.out)
