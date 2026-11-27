"""Exact pairwise certificates and realizable pure-state winners for a frozen pool."""
from itertools import product

import numpy as np

RT = 0.00198720425864083 * 298.15
SCORE_GAP_TOL = 0.001 + 1e-9  # two scores printed to 0.001 kcal/mol


def best_mask(a, b):
    """Occupied sites maximizing the weight ratio for a vs b (empty elsewhere)."""
    return a > b


def minimum_gap(ea, eb, a, b):
    """Infimum effective score of a minus b, when a has positive support."""
    take = best_mask(a, b)
    if np.any(b[take] == 0):
        return -np.inf
    return ea - eb - RT * (np.log(a[take]) - np.log(b[take])).sum()


def state_scores(energy, available, mask):
    with np.errstate(divide="ignore"):
        return energy - RT * np.log(available[:, mask]).sum(axis=1)


def assess(energy, available, rmsd):
    if (not np.isfinite(available).all() or np.any(available < 0)
            or np.any(available > 1) or not len(energy)):
        raise ValueError("Invalid archived pool")
    count, sites = available.shape
    excluded = np.zeros(count, dtype=bool)
    details = []
    masks = {tuple(np.zeros(sites, bool)), tuple(np.ones(sites, bool))}
    for site in range(sites):
        masks.add(tuple(np.arange(sites) == site))
    for i, a in enumerate(available):
        margins = np.full(count, -np.inf)
        for j, b in enumerate(available):
            if j != i:
                margins[j] = minimum_gap(energy[i], energy[j], a, b)
                masks.add(tuple(best_mask(a, b)))
        j = int(np.argmax(margins))
        excluded[i] = margins[j] > SCORE_GAP_TOL
        details.append((j if excluded[i] else None, margins[j]))
    if np.all(excluded):
        raise ArithmeticError("A finite dry winner must survive")
    lower = np.min(rmsd[~excluded])
    upper = np.inf
    winner_mask = None
    winners = None
    for mask in sorted(masks):
        scores = state_scores(energy, available, np.asarray(mask, dtype=bool))
        low = scores.min()
        if not np.isfinite(low):
            continue
        # Keep *all* candidates that could be selected under printed-score error.
        tied = np.flatnonzero(scores <= low + SCORE_GAP_TOL)
        error = rmsd[tied].max()
        if error < upper:
            upper, winner_mask, winners = error, mask, tied
    if lower > upper + 1e-10:
        raise ArithmeticError("Nonintersecting lower and feasible upper bounds")
    return details, lower, upper, winner_mask, winners


def check_enumeration():
    rng = np.random.default_rng(731)
    checks = 0
    for sites in range(7):
        for _ in range(8):
            a = rng.choice([0., 1e-90, .01, .4, .9, 1.], size=(4, sites))
            energy = rng.normal(size=4)
            levels = np.stack([state_scores(energy, a, np.asarray(mask, dtype=bool))
                               for mask in product((False, True), repeat=sites)])
            details, lower, upper, _, _ = assess(energy, a, np.arange(4.))
            for i in range(4):
                support = np.isfinite(levels[:, i])
                for j in range(4):
                    expected = (levels[support, i] - levels[support, j]).min()
                    if not np.isclose(minimum_gap(energy[i], energy[j], a[i], a[j]), expected):
                        raise ArithmeticError("Pair certificate differs from enumeration")
                    checks += 1
            for scores in levels:
                if np.isfinite(scores).any() and details[np.argmin(scores)][0] is not None:
                    raise AssertionError("Removed a positive-support winner")
    return checks
