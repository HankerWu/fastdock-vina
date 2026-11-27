"""Find cycles incompatible with a common free-ligand contribution to all rankings."""

import argparse
import csv
import itertools
import lzma
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from statistics import median

import numpy as np
from scipy.optimize import linprog

BENCH = Path(__file__).resolve().parents[2]
SOURCE = BENCH / "affinity_reference_ensemble/results"
PANELS = {"galectin": ("Gal1", "Gal3"), "kinase": ("AURKA", "CHEK1", "S6K1")}


def table(path):
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "rt") as stream:
        yield from csv.DictReader(stream)


def load_observations():
    kd, names = {}, {}
    for r in table(BENCH / "galectin_affinity/table1.csv"):
        i = int(r["compound"])
        for t in PANELS["galectin"]:
            kd[t, i] = 1000 * Fraction(r[f"{t}_Kd_uM"])
            names[i] = f"compound {i}"
    for r in table(BENCH / "torsion_kinase_affinity/panel.csv"):
        if r["target"] in PANELS["kinase"]:
            i = int(r["cid"])
            kd[r["target"], i] = Fraction(r["kd_nm"])
            names[i] = r["compound"]
    for r in table(SOURCE / "s6k_panel.csv"):
        if r["decision"] == "included":
            i = int(r["cid"])
            kd["S6K1", i] = Fraction(r["kd_nM"])
            names[i] = r["compound"]
    return kd, names


def score_data(kd):
    observations = defaultdict(list)
    unique = set()
    for name in ("primary_historical_cells.csv.xz", "s6k_cells.csv.xz"):
        for r in table(SOURCE / name):
            if int(r["nodes"]) != 1024:
                continue
            t, i, state = r["target"], int(r["ligand"]), r["state"]
            key = (t, i, state, r["seed"])
            if key in unique or kd[t, i] != Fraction(r["kd_nM"]):
                raise ValueError(("Duplicate or conflicting observation", key))
            unique.add(key)
            native = Fraction(r["native"])
            observations[t, i, "native", state].append(native)
            observations[t, i, "joint", state].append(native + Fraction(r["joint"]))
    states = defaultdict(list)
    for (t, i, arm, _), samples in observations.items():
        states[t, i, arm].append(median(samples))
    if len(unique) != 500:
        raise ValueError("Incomplete source")
    return {k: sum(s) / len(s) for k, s in states.items()}


def maximum_cycle_mean(edges):
    # A maximum-mean simple cycle visits at most n vertices. Dynamic programming
    # enumerates optimal fixed-length paths instead of listing all simple cycles.
    vertices = sorted({v for i, j, _, _ in edges for v in (i, j)})
    best, walk = None, []
    for start in vertices:
        paths = {start: (0, [])}
        for length in range(1, len(vertices) + 1):
            following = {}
            for edge in edges:
                i, j, _, weight = edge
                if i not in paths:
                    continue
                val = paths[i][0] + weight
                if j not in following or following[j][0] < val:
                    following[j] = val, paths[i][1] + [edge]
            paths = following
            if start in paths:
                val = paths[start][0] / length
                if best is None or val > best or (val == best and length < len(walk)):
                    best, walk = val, paths[start][1]
    if walk and len({e[0] for e in walk}) != len(walk):
        raise ArithmeticError("A shortest maximizer must be simple")
    return best, walk


def certify(edges, mean):
    vertices = sorted({v for i, j, _, _ in edges for v in (i, j)})
    col = {i: k for k, i in enumerate(vertices)}
    matrix = np.zeros((len(edges), len(vertices) + 1))
    rhs = []
    for k, (i, j, _, w) in enumerate(edges):
        matrix[k, col[i]], matrix[k, col[j]], matrix[k, -1] = 1, -1, -1
        rhs.append(-float(w))
    sol = linprog(
        [0.0] * len(vertices) + [1.0],
        A_ub=matrix,
        b_ub=rhs,
        bounds=[(0, 0)] + [(None, None)] * len(vertices),
        method="highs",
    )
    if mean is None:
        if sol.status != 3:
            raise ArithmeticError("Acyclic constraints allow arbitrary margins")
        return None, {}
    if not sol.success or abs(sol.fun - float(mean)) > 1e-8:
        raise ArithmeticError(("LP mismatch", sol.message, mean))
    # Exact dual certificate: the adjusted longest-path labels dominate every
    # edge, so no cycle in the full multigraph can have mean larger than mean.
    b = dict.fromkeys(vertices, 0)
    for _ in range(len(vertices) - 1):
        for i, j, _, w in edges:
            b[j] = max(b[j], b[i] + w - mean)
    if any(b[j] - b[i] < w - mean for i, j, _, w in edges):
        raise ArithmeticError("Primal potentials fail exact rational test")
    return sol.fun, b


def compatible_triangles(edges, removed):
    """Return strict conflicts on three ligands with no individually impossible unordered pair."""
    best_edge = {}
    for e in edges:
        key = e[:2]
        if key not in best_edge or e[3] > best_edge[key][3]:
            best_edge[key] = e
    vertices = sorted({v for pair in best_edge for v in pair})
    cycles = []
    for i, j, k in itertools.combinations(vertices, 3):
        if any(pair in removed for pair in ((i, j), (j, k), (i, k))):
            continue
        for triple in ((i, j, k, i), (i, k, j, i)):
            keys = list(itertools.pairwise(triple))
            if all(key in best_edge for key in keys):
                cycle = [best_edge[key] for key in keys]
                if sum(e[3] for e in cycle) >= 0:
                    cycles.append(cycle)
    return cycles


def output(path, records):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=records[0].keys(), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)


def run(dest):
    dest.mkdir(parents=True, exist_ok=True)
    kd, names = load_observations()
    scores = score_data(kd)
    summary, comparisons, witnesses, potentials = [], [], [], []
    for panel, targets in PANELS.items():
        for arm in ("native", "joint"):
            edges = []
            for t in targets:
                available = sorted(
                    i for target, i, model in scores if target == t and model == arm
                )
                for i, j in itertools.permutations(available, 2):
                    if kd[t, j] >= 3 * kd[t, i]:
                        edges.append((i, j, t, scores[t, i, arm] - scores[t, j, arm]))
            oriented = defaultdict(list)
            for i, j, _, w in edges:
                oriented[i, j].append(w)
            removed = {
                (i, j)
                for i, j in list(oriented)
                if (j, i) in oriented and max(oriented[i, j]) + max(oriented[j, i]) >= 0
            }
            for i, j, t, w in edges:
                comparisons.append(
                    {
                        "panel": panel,
                        "arm": arm,
                        "stronger": i,
                        "weaker": j,
                        "target": t,
                        "score_difference": str(w),
                        "removed": (i, j) in removed,
                    }
                )
            for mode in ("all", "without_impossible_reversals"):
                use = [e for e in edges if mode == "all" or e[:2] not in removed]
                mean, cycle = maximum_cycle_mean(use)
                lp, offsets = certify(use, mean)
                for i, b in sorted(offsets.items()):
                    potentials.append(
                        {
                            "panel": panel,
                            "arm": arm,
                            "graph": mode,
                            "ligand": i,
                            "common_offset": float(b),
                            "exact_offset": str(b),
                        }
                    )
                base = {"panel": panel, "arm": arm, "graph": mode}
                summary.append(
                    dict(
                        **base,
                        ligands=len({v for e in use for v in e[:2]}),
                        edges=len(use),
                        dropped_unordered_pairs=0
                        if mode == "all"
                        else len(removed) // 2,
                        max_mean=None if mean is None else float(mean),
                        exact_mean=None if mean is None else str(mean),
                        lp_max_mean=lp,
                        witness_length=len(cycle),
                        strictly_feasible=mean is None or mean < 0,
                    )
                )
                named = [("max_mean", cycle)]
                if mode == "without_impossible_reversals":
                    triangles = compatible_triangles(use, removed)
                    named.extend(
                        (f"compatible_triangle_{idx}", v)
                        for idx, v in enumerate(triangles, 1)
                    )
                for label, ordered in named:
                    for k, (i, j, t, w) in enumerate(ordered, 1):
                        witnesses.append(
                            dict(
                                **base,
                                certificate=label,
                                step=k,
                                target=t,
                                stronger=i,
                                weaker=j,
                                stronger_name=names[i],
                                weaker_name=names[j],
                                kd_stronger_nM=float(kd[t, i]),
                                kd_weaker_nM=float(kd[t, j]),
                                score_stronger=float(scores[t, i, arm]),
                                score_weaker=float(scores[t, j, arm]),
                                weight=float(w),
                                exact_weight=str(w),
                            )
                        )
    output(dest / "optimal_offsets.csv", potentials)
    output(dest / "summary.csv", summary)
    output(dest / "constraints.csv", comparisons)
    output(dest / "witnesses.csv", witnesses)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    run(p.parse_args().output)
