"""Joint feasibility of reference shifts and shared apo mixtures on frozen data."""

import argparse
import itertools
import math
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from statistics import median

import numpy as np
from scipy.optimize import linprog
from scipy.stats import kendalltau, spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "common_reference"))
from cycles import PANELS, SOURCE, load_observations, output, table

RT = 0.00198720425864083 * 298.15
TOL = 1e-8


def state_scores(kd, nodes):
    label = "primary" if nodes == 1024 else "sensitivity"
    s6k = "s6k_cells" if nodes == 1024 else "s6k_sensitivity"
    samples = defaultdict(list)
    seen = set()
    for path in [SOURCE / f"{label}_historical_cells.csv.xz", SOURCE / f"{s6k}.csv.xz"]:
        for row in table(path):
            if int(row["nodes"]) != nodes:
                continue
            t, i, r = row["target"], int(row["ligand"]), row["state"]
            key = t, i, r, row["seed"]
            if key in seen or kd[t, i] != Fraction(row["kd_nM"]):
                raise ValueError(f"Inconsistent input {key}")
            seen.add(key)
            n = Fraction(row["native"])
            samples["native", t, i, r].append(n)
            samples["joint", t, i, r].append(n + Fraction(row["joint"]))
    if len(seen) != 500:
        raise ValueError("Expected all 500 frozen poses")
    return {k: median(v) for k, v in samples.items()}


class Panel:
    def __init__(self, data, kd, targets, arm):
        self.targets = targets
        self.raw = {}
        edges = []
        self.kd = kd
        self.states = {}
        for t in targets:
            ids = sorted({i for a, target, i, r in data if a == arm and target == t})
            states = sorted({r for a, target, i, r in data if a == arm and target == t})
            self.states[t] = states
            for i in ids:
                self.raw[t, i] = [float(data[arm, t, i, r]) for r in states]
            edges.extend(
                (t, i, j)
                for i, j in itertools.permutations(ids, 2)
                if kd[t, j] >= 3 * kd[t, i]
            )
        self.edges = edges
        self.axes = [t for t in targets if len(self.states[t]) == 2]
        if len(self.axes) != 2 or any(
            len(s) not in (1, 2) for s in self.states.values()
        ):
            raise ValueError(self.states)

    def scores(self, point):
        p = dict(zip(self.axes, point, strict=True))
        values = {}
        for (t, i), arr in self.raw.items():
            if point == (None, None):
                v = sum(arr) / len(arr)
            elif len(arr) == 1 or p[t] == 0:
                v = arr[0]
            elif p[t] == 1:
                v = arr[1]
            else:
                a, b = arr
                v = a - RT * math.log1p(p[t] * math.expm1((a - b) / RT))
            values[t, i] = v
        return values

    def weights(self, point, edges):
        v = self.scores(point)
        return np.array([v[t, i] - v[t, j] for t, i, j in edges])


def solve(edges, w):
    ids = sorted({k for _, i, j in edges for k in (i, j)})
    col = {i: k for k, i in enumerate(ids)}
    a = np.zeros((len(edges), len(ids) + 1))
    for h, (_, i, j) in enumerate(edges):
        a[h, col[i]], a[h, col[j]], a[h, -1] = 1, -1, -1
    opt = linprog(
        [0.0] * len(ids) + [1.0],
        A_ub=a,
        b_ub=-w,
        bounds=[(0, 0)] + [(None, None)] * len(ids),
        method="highs",
    )
    if not opt.success:
        raise ArithmeticError(opt.message)
    return opt.fun, dict(zip(ids, opt.x[:-1], strict=True))


def positive_cycle(edges, weights):
    """Bellman-Ford on edge gains, followed by an explicit positivity check."""
    ids = {i for _, i, j in edges} | {j for _, i, j in edges}
    value = dict.fromkeys(ids, 0.0)
    pred = {}
    for _ in ids:
        recent = None
        for h, (_, i, j) in enumerate(edges):
            x = value[i] + weights[h]
            if x > value[j] + 1e-11:
                value[j], pred[j], recent = x, h, j
        if recent is None:
            return []
    for _ in ids:
        recent = edges[pred[recent]][1]
    tail, cycle = recent, []
    while True:
        h = pred[tail]
        cycle.append(h)
        tail = edges[h][1]
        if tail == recent:
            break
    cycle.reverse()
    if sum(weights[h] for h in cycle) <= TOL:
        return []
    return cycle


def box_values(panel, rect, edges):
    lo1, hi1, lo2, hi2 = rect
    # An edge belongs to just one target, hence only two opposing corners suffice.
    ends = [panel.weights((lo1, lo2), edges), panel.weights((hi1, hi2), edges)]
    return np.minimum(*ends)


def reversals(panel, weights):
    directed = defaultdict(list)
    for h, (_, i, j) in enumerate(panel.edges):
        directed[i, j].append(h)
    for (i, j), left in directed.items():
        if i >= j:
            continue
        for a, b in itertools.product(left, directed.get((j, i), [])):
            yield a, b, float(weights[a] + weights[b])


def population_only(panel, target, pairs):
    """With fixed offsets (all zero), ranking can change only at affine crossings."""
    factors = {}
    for (t, i), arr in panel.raw.items():
        if t == target:
            factors[i] = np.exp(-np.array([arr[0], arr[-1]]) / RT)
    delta = np.array([factors[i] - factors[j] for i, j in pairs])
    knots = {0.0, 1.0}
    for a, b in delta:
        if a * b < 0:
            knots.add(a / (a - b))
    # Each crossing tie is the average of the adjacent orderings. Endpoints
    # can tie without crossing and remain candidates in their own right.
    probes = [0.0, 1.0, *((a + b) / 2 for a, b in itertools.pairwise(sorted(knots)))]
    return max(float(((1 + np.sign(delta @ [1 - p, p])) / 2).sum()) for p in probes)


def classify(panel, edges):
    root = (0.0, 1.0, 0.0, 1.0)
    boxes, pending = [], [(root, 0)]
    while pending:
        rect, depth = pending.pop()
        lower = box_values(panel, rect, edges)
        cycle = positive_cycle(edges, lower)
        if cycle:
            boxes.append((rect, cycle, lower))
            continue
        x, X, y, Y = rect
        points = [((x + X) / 2, (y + Y) / 2), *itertools.product((x, X), (y, Y))]
        for point in points:
            w = panel.weights(point, edges)
            value, offsets = solve(edges, w)
            if value < -TOL:
                return "feasible", [], (point, offsets, value)
        if depth == 22:
            raise ArithmeticError(
                "Reached the declared box precision without a certificate"
            )
        if X - x >= Y - y:
            mid = (x + X) / 2
            children = [(x, mid, y, Y), (mid, X, y, Y)]
        else:
            mid = (y + Y) / 2
            children = [(x, X, y, mid), (x, X, mid, Y)]
        pending.extend((r, depth + 1) for r in reversed(children))
    return "infeasible", boxes, None


def run(dest):
    dest.mkdir(parents=True, exist_ok=True)
    kd, names = load_observations()
    results = defaultdict(list)
    for nodes in (1024, 8192):
        data = state_scores(kd, nodes)
        for family, targets in PANELS.items():
            for arm in ("native", "joint"):
                panel = Panel(data, kd, targets, arm)
                base = {"nodes": nodes, "family": family, "arm": arm}
                for (t, i), vals in panel.raw.items():
                    for r, v in zip(panel.states[t], vals, strict=True):
                        results["states"].append(
                            dict(
                                **base,
                                target=t,
                                ligand=i,
                                name=names[i],
                                kd_nM=float(kd[t, i]),
                                state=r,
                                energy=v,
                            )
                        )
                equal = panel.scores((0.5, 0.5))
                for t in targets:
                    ids = sorted(i for target, i in equal if target == t)
                    observed = [float(kd[t, i]) for i in ids]
                    predicted = [equal[t, i] for i in ids]
                    pairs = [(i, j) for target, i, j in panel.edges if t == target]
                    credit = sum(
                        (1 + np.sign(equal[t, j] - equal[t, i])) / 2 for i, j in pairs
                    )
                    results["baseline_ranks"].append(
                        dict(
                            **base,
                            target=t,
                            ligands=len(ids),
                            pairs=len(pairs),
                            spearman=spearmanr(observed, predicted).statistic,
                            kendall_tau_b=kendalltau(observed, predicted).statistic,
                            correct=credit,
                            population_only_max_credit=population_only(panel, t, pairs),
                        )
                    )
                removed = set()
                lower = box_values(panel, (0, 1, 0, 1), panel.edges)
                rev = list(reversals(panel, lower))
                if any(abs(val) < TOL for a, b, val in rev):
                    raise ArithmeticError("Boundary reversal requires tie handling")
                for a, b, val in rev:
                    t, i, j = panel.edges[a]
                    u, _, _ = panel.edges[b]
                    if val > 0:
                        removed.add((i, j))
                    results["reversals"].append(
                        dict(
                            **base,
                            i=i,
                            j=j,
                            target_i_stronger=t,
                            target_j_stronger=u,
                            min_cycle_sum=val,
                            always_conflicts=int(val > 0),
                        )
                    )
                remaining = [
                    e for e in panel.edges if tuple(sorted(e[1:])) not in removed
                ]
                for graph, edges in [
                    ("all", panel.edges),
                    ("without_unavoidable_pairs", remaining),
                ]:
                    prefix = dict(**base, graph=graph)
                    status, boxes, feasible = classify(panel, edges)
                    uniform, _ = solve(edges, panel.weights((0.5, 0.5), edges))
                    arithmetic, _ = solve(edges, panel.weights((None, None), edges))
                    floor = None
                    for k, (rect, cycle, lo) in enumerate(boxes):
                        strength = sum(lo[h] for h in cycle) / len(cycle)
                        floor = strength if floor is None else min(floor, strength)
                        for order, h in enumerate(cycle):
                            t, i, j = edges[h]
                            results["cycles"].append(
                                dict(
                                    **prefix,
                                    box=k,
                                    xlo=rect[0],
                                    xhi=rect[1],
                                    ylo=rect[2],
                                    yhi=rect[3],
                                    step=order,
                                    target=t,
                                    stronger=i,
                                    weaker=j,
                                    edge_lower=float(lo[h]),
                                )
                            )
                    if feasible:
                        point, offsets, _ = feasible
                        for i, b in offsets.items():
                            results["feasible"].append(
                                dict(
                                    **prefix, x=point[0], y=point[1], ligand=i, offset=b
                                )
                            )
                    results["summary"].append(
                        dict(
                            **prefix,
                            first_axis=panel.axes[0],
                            second_axis=panel.axes[1],
                            edges=len(edges),
                            reverse_comparisons=len(rev),
                            removed_unordered_pairs=len(removed),
                            equal_energy_mean_minimax=arithmetic,
                            equal_population_minimax=uniform,
                            status=status,
                            proof_boxes=len(boxes),
                            certified_lower_bound=floor,
                            witness_upper_bound=feasible[2] if feasible else None,
                        )
                    )
    for name, rows in results.items():
        output(dest / f"{name}.csv", rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dest", type=Path)
    run(parser.parse_args().dest)
