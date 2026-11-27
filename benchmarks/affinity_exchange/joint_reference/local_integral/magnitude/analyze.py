"""Global additive-reference amplitude bounds for two- and three-target panels."""
import argparse
import csv
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

HERE = Path(__file__).resolve().parent
RT = 0.00198720425864083 * 298.15
TOL = 0.0009


def inputs(nodes, family, arm, observation):
    parent = HERE.parent
    states = defaultdict(dict)
    kd = {}
    with (parent / "results/states.csv").open() as stream:
        for r in csv.DictReader(stream):
            if (r["nodes"], r["family"], r["arm"]) != (
                str(nodes), family, "native" if arm == "zero" else arm
            ):
                continue
            key = r["target"], int(r["ligand"])
            states[key][r["state"]] = 0. if arm == "zero" else float(r["energy"])
            kd[key] = float(r["kd_nM"])
    sem = defaultdict(float)
    if observation == "one_sem":
        with (HERE.parents[3] / "galectin_affinity/table1.csv").open() as stream:
            for r in csv.DictReader(stream):
                for t in ("Gal1", "Gal3"):
                    key = t, int(r["compound"])
                    assert kd[key] == 1000 * float(r[t + "_Kd_uM"])
                    sem[key] = 1000 * float(r[t + "_error_uM"])
    keys = sorted(states)
    targets = sorted({t for t, i in keys})
    structures = {t: sorted(set.union(*(set(states[key]) for key in keys if key[0] == t)))
                  for t in targets}
    energies = []
    for key in keys:
        ids = structures[key[0]]
        assert set(states[key]) == set(ids) and 1 <= len(ids) <= 2
        energies.append([states[key][ids[0]], states[key][ids[-1]]])
    lo = RT * np.log([1e-9 * (kd[k] - sem[k]) for k in keys])
    hi = RT * np.log([1e-9 * (kd[k] + sem[k]) for k in keys])
    return keys, structures, np.array(energies), lo, hi


def cycles(keys, targets):
    """Each row is one undirected simple cycle. Closing the walk cancels all offsets."""
    ligands = {t: {i for u, i in keys if u == t} for t in targets}
    index = {key: n for n, key in enumerate(keys)}
    named = []
    for size in range(2, len(targets) + 1):
        for ts in itertools.combinations(targets, size):
            if size == 2:
                choices = itertools.combinations(sorted(ligands[ts[0]] & ligands[ts[1]]), 2)
            else:
                assert size == 3
                choices = itertools.product(*(sorted(ligands[t] & ligands[ts[(j+1) % size]])
                                              for j, t in enumerate(ts)))
            for lig in choices:
                if len(set(lig)) == size:
                    named.append([(t, lig[j], lig[j-1]) for j, t in enumerate(ts)])
    # The extra zero column represents absent targets in cycles on target subsets.
    plus = np.full((len(named), len(targets)), len(keys), dtype=int)
    minus = plus.copy()
    for n, cycle in enumerate(named):
        for t, i, j in cycle:
            k = targets.index(t)
            plus[n, k], minus[n, k] = index[t, i], index[t, j]
    return named, plus, minus


def solve(nodes, family, arm, observation):
    keys, structures, energy, obs_lo, obs_hi = inputs(nodes, family, arm, observation)
    targets = list(structures)
    axes = [t for t in targets if len(structures[t]) == 2]
    axis = [axes.index(t) if t in axes else len(axes) for t, _ in keys]
    named, plus, minus = cycles(keys, targets)
    length = np.array([2*len(c) for c in named])
    def diff(v1, v2):
        return np.append(v1, 0)[plus] - np.append(v2, 0)[minus]
    obs_min = diff(obs_lo, obs_hi).sum(axis=1)
    obs_max = diff(obs_hi, obs_lo).sum(axis=1)

    def mixed(p):
        # log1p/expm1 preserve the pure-state limit and resolve close energies.
        pp = np.append(p, 0)[axis]
        return energy[:, 0] - RT * np.log1p(pp * np.expm1((energy[:, 0] - energy[:, 1]) / RT))

    def bound(left, right):
        el, er = mixed(left), mixed(right)
        dl, dr = diff(el, el), diff(er, er)
        low = (obs_min - np.maximum(dl, dr).sum(axis=1)) / length
        neg_high = (np.minimum(dl, dr).sum(axis=1) - obs_max) / length
        r = np.maximum(low, neg_high)
        n = int(np.argmax(r))
        return max(0., float(r[n])), named[n]

    # Use feasible points to improve the upper bound. Every terminal box still
    # needs its own valid lower bound, regardless of point selection.
    points = list(itertools.product((0., .5, 1.), repeat=len(axes)))
    best, p = min((bound(v, v)[0], v) for v in points)
    p = np.array(p)
    boxes = []

    def visit(left, right, path):
        nonlocal best, p
        lower, cycle = bound(left, right)
        center = (left + right) / 2
        value, _ = bound(center, center)
        if value < best:
            best, p = value, center
        if lower >= best - TOL:
            boxes.append(dict(path=path, lower=math.floor((lower - 1e-10) * 1e6) / 1e6, cycle=cycle))
            return
        dim = len(path) % len(axes)
        middle = center[dim]
        new_right = right.copy()
        new_right[dim] = middle
        new_left = left.copy()
        new_left[dim] = middle
        visit(left, new_right, path+"0")
        visit(new_left, right, path+"1")

    full = (np.zeros(len(axes)), np.ones(len(axes)))
    root_bound, root_cycle = bound(*full)
    full_domain = dict(path="", lower=math.floor((root_bound - 1e-10) * 1e6) / 1e6, cycle=root_cycle)
    if family == "kinase" and arm == "joint":
        # Post-hoc single-cycle check. A common positive scale preserves the
        # endpoint signs, even if applied before taking the receptor mixture.
        obs = lo_score = hi_score = 0.
        for t, i, j in root_cycle:
            a, b = keys.index((t, i)), keys.index((t, j))
            assert obs_lo[a] == obs_hi[a] and obs_lo[b] == obs_hi[b]
            obs += obs_lo[a] - obs_hi[b]
            ends = energy[a] - energy[b]
            lo_score += min(ends)
            hi_score += max(ends)
        if obs < 0. < lo_score or hi_score < 0. < obs:
            full_domain["positive_scale_lower"] = math.floor(abs(obs)/(2*len(root_cycle))*1e6)/1e6
    visit(*full, "")
    ligands = sorted({i for _, i in keys})
    design = np.zeros((len(keys), len(targets) + len(ligands)))
    for n, (t, i) in enumerate(keys):
        design[n, targets.index(t)] = design[n, len(targets) + ligands.index(i)] = 1.
    s = mixed(p)
    result = linprog([0.]*design.shape[1]+[1.],
                     A_ub=np.column_stack((np.vstack((design, -design)), -np.ones(2*len(keys)))),
                     b_ub=np.r_[obs_hi-s, s-obs_lo],
                     bounds=[(0., 0.)]+[(None,None)]*(design.shape[1]-1)+[(0.,None)],
                     method="highs")
    if not result.success or abs(result.fun-best) > 1e-7:
        raise ValueError((result.message, best, result.fun))
    deviations = np.maximum(0., np.maximum(s + design @ result.x[:-1] - obs_hi,
                                            obs_lo - s - design @ result.x[:-1]))
    upper = math.ceil((max(deviations) + 1e-10) * 1e6) / 1e6
    lower = min(b["lower"] for b in boxes)
    if upper-lower > .001:
        raise ArithmeticError((upper, lower))
    return dict(nodes=nodes, family=family, arm=arm, observation=observation,
                cells=len(keys), cycles=len(named), lower=lower, upper=upper,
                structures=structures, axes=axes, population=p.tolist(),
                target_offsets=dict(zip(targets, result.x[:len(targets)], strict=True)),
                ligand_offsets=dict(zip(ligands, result.x[len(targets):-1], strict=True)),
                full_domain=full_domain, boxes=boxes)


def main(dest):
    dest.mkdir(parents=True, exist_ok=True)
    reports = []
    for nodes in (1024, 8192):
        for family in ("galectin", "kinase"):
            arms = ("zero", "native", "joint") if nodes == 1024 else ("joint",)
            for arm in arms:
                for obs in (("mean", "one_sem") if family == "galectin" else ("mean",)):
                    reports.append(solve(nodes, family, arm, obs))
    (dest / "witnesses.json").write_text(json.dumps(reports, indent=1) + "\n")
    with (dest / "summary.csv").open("w", newline="") as stream:
        fields = ("nodes", "family", "arm", "observation", "cells", "cycles", "lower", "upper")
        out = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        out.writeheader()
        out.writerows(reports)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    main(parser.parse_args().output)
