"""Exact state enumeration and primal/dual certificates for the open pool bounds."""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.special import logsumexp

from envelope import RT
from study import KEY

TOL = 1e-8
ROUNDING = .001


def columns(energies, avail):
    n = avail.shape[1]
    states = np.arange(2**n, dtype=np.uint32)
    occupied = (states[:, None] >> np.arange(n)) & 1
    # An occupied zero-availability site kills the kernel. Avoid 0 * log(0).
    logs = -energies[:, None] / RT + np.log(np.where(avail > 0, avail, 1)) @ occupied.T
    logs[((avail == 0).astype(int) @ occupied.T) > 0] = -np.inf
    return logs


def game(logs, target, favor):
    """Normalize each nondead column; reconstruct the original full-site law."""
    weight = logs.copy()
    weight[target] += favor * ROUNDING / RT
    peak = weight.max(axis=0)
    live = np.isfinite(peak)
    valid = np.flatnonzero(live)
    weights = np.exp(weight[:, live] - peak[live])
    rivals = np.delete(np.arange(len(weights)), target)
    margins = weights[target] - weights[rivals]
    count = len(rivals)
    # Dual game: choose one mixture of opponents that bounds every configuration.
    result = linprog(np.r_[np.zeros(count), 1.],
                     A_ub=np.column_stack((margins.T, -np.ones(len(valid)))),
                     b_ub=np.zeros(len(valid)), A_eq=[np.r_[np.ones(count), 0.]],
                     b_eq=[1.], bounds=[(0, None)] * count + [(None, None)],
                     method="highs", options={"primal_feasibility_tolerance": 1e-10,
                                             "dual_feasibility_tolerance": 1e-10})
    if not result.success:
        raise ValueError(result.message)
    dual = np.maximum(result.x[:-1], 0)
    dual /= dual.sum()
    upper = (dual @ margins).max()
    primal = np.maximum(-result.ineqlin.marginals, 0)
    primal /= primal.sum()
    lower = (margins @ primal).min()
    support = np.flatnonzero(primal > 0)
    logq = np.log(primal[support]) - peak[valid[support]]
    logq -= logsumexp(logq)
    # Read back the real apo weights against original kernels, independently of LP scaling.
    adjusted = logsumexp(weight[:, valid[support]] + logq, axis=1)
    if lower > TOL and adjusted[target] <= max(adjusted[rivals]):
        raise ArithmeticError("Reconstructed water law changes winner")
    if lower > upper + TOL or upper - lower > 2e-7:
        raise ArithmeticError("Primal-dual gap at prediction boundary")
    return dict(lower=float(lower), upper=float(upper), q_states=valid[support].tolist(),
                q=np.exp(logq).tolist(), logq=logq.tolist(),
                rival_weights={str(int(j)+1):float(z) for j,z in zip(rivals,dual) if z>0})


def main(source, out):
    w, c = time.perf_counter(), time.process_time()
    cells = pd.read_csv(source / "cells.csv.gz", float_precision="round_trip")
    modes = pd.read_csv(source / "modes.csv.gz", float_precision="round_trip")
    factors = pd.read_csv(source / "factors.csv.gz", float_precision="round_trip")
    keys = ["order", *KEY]
    grouped_modes, grouped_factors = modes.groupby(keys), factors.groupby(keys)
    decisions, evidence, bounds = [], [], []
    for key, table in cells.groupby(keys):
        cell = table.iloc[0]
        if cell.closed:
            continue
        prefix = {k: v.item() if isinstance(v, np.generic) else v for k, v in zip(keys, key, strict=True)}
        m = grouped_modes.get_group(key).sort_values("model")
        inputs = grouped_factors.get_group(key).pivot(index="model", columns="site",values="availability")
        if not np.array_equal(m.model, inputs.index):
            raise ValueError("Mode index mismatch")
        logs = columns(m.dry.to_numpy(), inputs.to_numpy())
        lower = np.inf
        upper = cell.witness_upper
        for pos, row in enumerate(m.itertuples()):
            if row.dominated:
                continue
            best = game(logs, pos, +1)
            worst = game(logs, pos, -1)
            label = ("impossible" if best["upper"] < -TOL else
                     "realizable" if worst["lower"] > TOL else "unresolved")
            if label != "impossible":
                lower = min(lower, row.rmsd)
            if label == "realizable":
                upper = min(upper, row.rmsd)
            decisions.append(dict(**prefix, model=row.model, rmsd=row.rmsd,
                                  classification=label, favorable_dual_upper=best["upper"],
                                  adverse_primal_lower=worst["lower"]))
            evidence.append(dict(**prefix, model=row.model, sites=inputs.columns.to_list(),
                                 favorable=best, adverse=worst))
        if lower > upper + 1e-10:
            raise ArithmeticError("Inconsistent envelope")
        bounds.append(dict(**prefix, pair_lower=cell.best_possible_lower,
                           pure_upper=cell.witness_upper, lower=lower, upper=upper,
                           oracle=cell.oracle, closed=abs(upper-lower)<1e-10))
    out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(decisions).to_csv(out / "game_modes.csv",index=False)
    pd.DataFrame(bounds).to_csv(out / "game_pools.csv",index=False)
    (out / "game_witnesses.json").write_text(json.dumps(evidence,indent=1)+"\n")
    (out / "game_cost.json").write_text(json.dumps(dict(wall_s=time.perf_counter()-w,
                                                       cpu_s=time.process_time()-c),indent=2)+"\n")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("source",type=Path)
    p.add_argument("output",type=Path)
    args=p.parse_args()
    main(args.source,args.output)
