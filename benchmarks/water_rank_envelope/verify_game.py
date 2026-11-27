"""Independent direct-product enumeration of all original-law and rival certificates."""
import argparse
import collections
import csv
import gzip
import json
import time
from pathlib import Path

import numpy as np

KEYS = ("order", "panel", "target", "cap", "query", "receiver", "seed", "search")
RT = np.longdouble("0.00198720425864083") * np.longdouble("298.15")
ROUND = np.longdouble("0.001")
CUTOFF = 1e-8


def key(row):
    return tuple(str(row[name]) for name in KEYS)


def rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return list(csv.DictReader(f))


def direct_product(energy, availability):
    """Double the column table once per site; bit 0 corresponds to its first column."""
    weight = np.exp((energy.min() - energy)[:, None] / RT)
    for site in availability.T:
        occupied = weight * site[:, None]
        if np.any((occupied == 0) & (weight != 0) & (site[:, None] != 0)):
            raise FloatingPointError("Nonzero kernel underflowed")
        weight = np.concatenate((weight, occupied), axis=1)
    return weight


def check_law(w, target, entry):
    states = np.array(entry["q_states"])
    q = np.array(entry["q"], dtype=np.longdouble)
    assert len(states) == len(set(states)) == len(q)
    assert (states >= 0).all() and (states < w.shape[1]).all() and (q > 0).all()
    rivals = [i for i in range(len(w)) if i != target]
    alpha = np.array([entry["rival_weights"].get(str(i + 1), 0) for i in rivals], dtype=np.longdouble)
    assert (alpha >= 0).all() and set(map(int, entry["rival_weights"])) <= {i + 1 for i in rivals}
    law_error = max(abs(q.sum()-1), abs(alpha.sum()-1), np.max(abs(q-np.exp(entry["logq"]))))
    peak = w.max(axis=0)
    active = peak > 0
    # A common positive column factor cancels from a sign test. Reconstruct the
    # original lower margin *without* first normalizing the emitted water law.
    value = w[:, states] @ q
    primal = (value[target] - value[rivals].max()) / (peak[states] @ q)
    rival_mixture = alpha @ w[rivals]
    dual = ((w[target, active] - rival_mixture[active]) / peak[active]).max()
    residual = max(abs(primal-entry["lower"]), abs(dual-entry["upper"]))
    assert primal <= dual + 1e-10
    return float(primal), float(dual), float(law_error), float(residual)


def run(path):
    mode = collections.defaultdict(dict)
    factors = collections.defaultdict(dict)
    for r in rows(path/"modes.csv.gz"):
        mode[key(r)][int(r["model"])] = r
    for r in rows(path/"factors.csv.gz"):
        factors[key(r)][int(r["model"]), r["site"]] = np.longdouble(r["availability"])
    bounds = {key(r): r for r in rows(path/"game_pools.csv")}
    open_cells = {key(r) for r in rows(path/"cells.csv.gz") if r["closed"] == "False"}
    assert open_cells == bounds.keys()
    decisions = {(key(r), int(r["model"])): r for r in rows(path/"game_modes.csv")}
    witnesses = collections.defaultdict(list)
    for e in json.loads((path/"game_witnesses.json").read_text()):
        witnesses[key(e)].append(e)
    assert witnesses.keys() == bounds.keys()

    labels, checked = collections.Counter(), set()
    residual = lawerr = pure_residual = 0.
    state_count = model_state_count = 0
    per_pool = []
    for k, certificate in bounds.items():
        m = mode[k]
        ids = sorted(m)
        assert ids == list(range(1, len(m)+1))
        site = sorted({s for _, s in factors[k]})
        energy = np.array([m[i]["dry"] for i in ids], dtype=np.longdouble)
        rmsd = np.array([m[i]["rmsd"] for i in ids], dtype=np.longdouble)
        a = np.array([[factors[k][i, s] for s in site] for i in ids])
        w = direct_product(energy, a)
        live = w.max(axis=0) > 0
        weights = w[:, live]
        # Exact best pure-state worst winner, including all rounding-compatible ties.
        possible = weights >= weights.max(axis=0)[None] * np.exp(-(ROUND+1e-9)/RT)
        worst = np.where(possible, rmsd[:, None], -np.inf).max(axis=0)
        best_pure = float(worst.min())
        pure_residual = max(pure_residual, abs(best_pure-float(certificate["pure_upper"])))
        assert best_pure >= float(certificate["lower"]) - 1e-12
        assert best_pure <= float(certificate["pure_upper"]) + 1e-12
        state_count += w.shape[1]
        per_pool.append(dict(zip(KEYS, k)) | dict(sites=len(site), states=w.shape[1],
                        all_dead=int((~live).sum()), exact_best_pure_A=best_pure))
        for e in witnesses[k]:
            model = int(e["model"])
            assert (k, model) not in checked and e["sites"] == site
            checked.add((k, model))
            limits = {}
            for sign, name in ((1, "favorable"), (-1, "adverse")):
                shifted = w.copy()
                shifted[model-1] *= np.exp(sign*ROUND/RT)
                lo, hi, prob_err, margin_err = check_law(shifted, model-1, e[name])
                model_state_count += w.shape[1]
                residual, lawerr = max(residual, margin_err), max(lawerr, prob_err)
                limits[name] = lo, hi
            label = ("impossible" if limits["favorable"][1] < -CUTOFF else
                     "realizable" if limits["adverse"][0] > CUTOFF else "unresolved")
            row = decisions[k, model]
            assert label == row["classification"]
            assert abs(float(row["favorable_dual_upper"])-limits["favorable"][1]) < 1e-10
            assert abs(float(row["adverse_primal_lower"])-limits["adverse"][0]) < 1e-10
            labels[label] += 1
    assert checked == decisions.keys() and residual < 1e-10 and lawerr < 1e-10
    return dict(pools=len(bounds), certificates=len(checked), favorable_and_adverse_checks=2*len(checked),
                states_across_pools=state_count, states_across_certificate_signs=model_state_count,
                counts=labels, max_abs_probability_error=lawerr, max_abs_game_margin_error=residual,
                max_abs_exact_pure_minus_previous_witness_A=pure_residual, per_pool=per_pool,
                accumulator=str(np.finfo(np.longdouble).dtype))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    tic, cpu = time.perf_counter(), time.process_time()
    result = run(args.results)
    result.update(wall_s=time.perf_counter()-tic, cpu_s=time.process_time()-cpu)
    args.out.write_text(json.dumps(result, indent=2)+"\n")
