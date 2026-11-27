"""Verify published answers by direct association ratios and 50-digit arithmetic."""

import argparse
import csv
import itertools
import json
import lzma
from collections import defaultdict
from decimal import Decimal as D
from decimal import getcontext
from functools import partial
from pathlib import Path
from statistics import median

getcontext().prec = 50
RT = D("0.00198720425864083") * D("298.15")
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[2] / "affinity_reference_ensemble/results"


def read(path):
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "rt") as stream:
        return list(csv.DictReader(stream))


def read_raw(nodes):
    primary = nodes == "1024"
    files = (
        "primary_historical_cells" if primary else "sensitivity_historical_cells",
        "s6k_cells" if primary else "s6k_sensitivity",
    )
    energies = defaultdict(list)
    kds = {}
    unique = set()
    for name in files:
        for row in read(SOURCE / (name + ".csv.xz")):
            if row["nodes"] != nodes:
                continue
            key = (row["target"], int(row["ligand"]))
            kds[key] = D(row["kd_nM"])
            state = row["state"]
            unique.add((*key, state, row["seed"]))
            n = D(row["native"])
            energies["native", *key, state].append(n)
            energies["joint", *key, state].append(n + D(row["joint"]))
    assert len(unique) == 500
    aggr = {key: median(values) for key, values in energies.items()}
    return aggr, kds


def grouped(rows, fields):
    result = defaultdict(list)
    for row in rows:
        result[tuple(row[f] for f in fields)].append(row)
    return result


def compare(energies, arm, states, axes, edge, point):
    t, i, j = edge
    if len(states[t]) == 1:
        state = states[t][0]
        return energies[arm, t, i, state] - energies[arm, t, j, state]
    p = point[axes.index(t)]

    def strength(ligand):
        return sum(
            w * (-energies[arm, t, ligand, r] / RT).exp()
            for r, w in zip(states[t], (1 - p, p), strict=True)
        )

    return -RT * (strength(i) / strength(j)).ln()


def verify_native_sources(raw):
    """Link the 500-cell path to the independent original native readback route."""
    parent = HERE.parent
    order = json.loads((parent / "result.json").read_text())["states"]
    entries = read(parent / "medians.csv")
    assert len(entries) == 73
    for energies, kd in raw.values():
        expected = set()
        for row in entries:
            target, ligand = row["target"], int(row["ligand"])
            assert kd[target, ligand] == D(row["kd_nM"])
            for k, state in enumerate(order[target], start=1):
                key = ("native", target, ligand, state)
                expected.add(key)
                assert energies[key] == D(row[f"score_{k}"])
        assert expected == {key for key in energies if key[0] == "native"}


def check(dest):
    fields = ("nodes", "family", "arm", "graph")
    cycles = grouped(read(dest / "cycles.csv"), fields)
    fits = grouped(read(dest / "feasible.csv"), fields)
    revs = grouped(read(dest / "reversals.csv"), fields[:3])
    raw = {n: read_raw(n) for n in ("1024", "8192")}
    verify_native_sources(raw)
    for summary in read(dest / "summary.csv"):
        key = tuple(summary[f] for f in fields)
        nodes, family, arm, graph = key
        energies, kd = raw[nodes]
        axes = (summary["first_axis"], summary["second_axis"])
        targets = (*axes, "S6K1") if family == "kinase" else axes
        states = {
            t: sorted({r for model, tgt, _, r in energies if tgt == t}) for t in targets
        }
        edges = set()
        for t in targets:
            ligands = [i for tgt, i in kd if tgt == t]
            edges.update(
                (t, i, j)
                for i, j in itertools.permutations(ligands, 2)
                if kd[t, j] >= 3 * kd[t, i]
            )

        gap = partial(compare, energies, arm, states, axes)
        removed = set()
        for row in revs[key[:3]]:
            i, j = int(row["i"]), int(row["j"])
            a, b = (row["target_i_stronger"], i, j), (row["target_j_stronger"], j, i)
            assert a in edges and b in edges
            # Enumerate the four independent pure choices here, unlike the
            # main code, which bounds each single-target ratio separately.
            lowest = min(
                gap(a, xy) + gap(b, xy)
                for xy in itertools.product([D(0), D(1)], repeat=2)
            )
            assert abs(lowest - D(row["min_cycle_sum"])) < D("1e-10")
            assert (lowest > 0) == bool(int(row["always_conflicts"]))
            if lowest > 0:
                removed.add(tuple(sorted((i, j))))
        assert len(removed) == int(summary["removed_unordered_pairs"])
        if graph != "all":
            edges = {e for e in edges if tuple(sorted(e[1:])) not in removed}
        assert len(edges) == int(summary["edges"])

        if summary["status"] == "feasible":
            rows = fits[key]
            point = (D(rows[0]["x"]), D(rows[0]["y"]))
            offsets = {int(row["ligand"]): D(row["offset"]) for row in rows}
            assert {i for _, i, j in edges} | {j for _, i, j in edges} == offsets.keys()
            violation = max(
                gap(e, point) + offsets[e[1]] - offsets[e[2]] for e in edges
            )
            assert violation < 0 and abs(
                violation - D(summary["witness_upper_bound"])
            ) < D("1e-10")
            # Strict feasibility also holds just inside pure-state endpoints.
            interior = tuple(min(D(".999"), max(D(".001"), v)) for v in point)
            interior_margin = min(
                -gap(e, interior) - offsets[e[1]] + offsets[e[2]] for e in edges
            )
            assert interior_margin > 0
            print(
                nodes,
                family,
                arm,
                graph,
                "min interior margin",
                round(interior_margin, 8),
            )
            continue

        boxes = []
        strengths = []
        for items in grouped(cycles[key], ["box"]).values():
            first = items[0]
            bounds = tuple(D(first[k]) for k in ["xlo", "xhi", "ylo", "yhi"])
            x, X, y, Y = bounds
            boxes.append(bounds)
            assert 0 <= x < X <= 1 and 0 <= y < Y <= 1
            points = list(itertools.product((x, X), (y, Y)))
            steps = [(r["target"], int(r["stronger"]), int(r["weaker"])) for r in items]
            assert set(steps) <= edges
            assert all(a[2] == b[1] for a, b in zip(steps, steps[1:] + steps[:1]))
            bounds_sum = D(0)
            for e, r in zip(steps, items, strict=True):
                value = min(gap(e, p) for p in points)
                assert abs(value - D(r["edge_lower"])) < D("1e-10")
                bounds_sum += value
            assert bounds_sum > 0
            strengths.append(bounds_sum / len(steps))
        assert abs(min(strengths) - D(summary["certified_lower_bound"])) < D("1e-10")
        assert sum((X - x) * (Y - y) for x, X, y, Y in boxes) == 1
        for a, b in itertools.combinations(boxes, 2):
            assert a[0] >= b[1] or b[0] >= a[1] or a[2] >= b[3] or b[2] >= a[3]
    print("All reversal, rectangle, cycle and direct ordering certificates verified.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    check(parser.parse_args().results)
