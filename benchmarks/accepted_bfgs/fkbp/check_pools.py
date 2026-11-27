"""Read AD4 coordinates and scores without using the search or scoring analysis code."""

import argparse
import csv
import hashlib
import itertools
import json
import math
import statistics
import tarfile
from decimal import Decimal
from pathlib import Path

RECEPTORS = ("1FKG", "1FKH")
SEEDS = (131071, 196613, 262147)


def ranks(xs):
    return [1 + sum(v < x for v in xs) + (sum(v == x for v in xs) - 1) / 2 for x in xs]


def kendall(xs, ys):
    sx = sy = net = 0
    for i, j in itertools.combinations(range(len(xs)), 2):
        a = (xs[i] > xs[j]) - (xs[i] < xs[j])
        b = (ys[i] > ys[j]) - (ys[i] < ys[j])
        net += a * b
        sx += a * a
        sy += b * b
    return net / math.sqrt(sx * sy)


def pool(data):
    score = None
    models = []
    for block in data.split(b"ENDMDL"):
        lines = block.splitlines()
        rows = [l for l in lines if l.startswith((b"ATOM  ", b"HETATM"))]
        if not rows:
            continue
        if score is None:
            score = Decimal(
                next(
                    l.split()[3].decode()
                    for l in lines
                    if l.startswith(b"REMARK VINA RESULT:")
                )
            )
        entries = sorted((int(l[6:11]), l[30:54].strip()) for l in rows)
        assert len({i for i, c in entries}) == len(entries)
        h = hashlib.sha256(
            b"\n".join(str(i).encode() + b" " + coord for i, coord in entries)
        ).hexdigest()
        models.append(h)
    assert score is not None and models
    return score, models


def check(work, benchmark):
    with (benchmark / "human_panel/affinity/ligands.csv").open() as f:
        bounds = {
            int(r["compound"]): (Decimal(r["ki_low_nm"]), Decimal(r["ki_high_nm"]))
            for r in csv.DictReader(f)
        }
    edges = [
        (a, b)
        for a, b in itertools.permutations(bounds, 2)
        if bounds[b][0] >= 3 * bounds[a][1]
    ]
    point_ids = sorted(k for k, v in bounds.items() if v[0] == v[1])
    assert len(edges) == 15 and len(point_ids) == 7
    arms = ("official", "parent", "fixed")
    scores = {}
    coords = {}
    with tarfile.open(benchmark / "fkbp_ad4/pools.tar.xz") as official:
        for comp, rec, seed in itertools.product(sorted(bounds), RECEPTORS, SEEDS):
            name = f"pools/{comp}_{rec}_{seed}_ad4.pdbqt"
            for arm in arms:
                data = (
                    official.extractfile(name).read()
                    if arm == "official"
                    else (work / arm / name).read_bytes()
                )
                scores[arm, comp, rec, seed], coords[arm, comp, rec, seed] = pool(data)
    results = {}
    flags = {}
    for arm in arms:
        receptor = {
            r: {
                c: statistics.median(scores[arm, c, r, s] for s in SEEDS)
                for c in bounds
            }
            for r in RECEPTORS
        }
        receptor["receptor_mean"] = {
            c: statistics.mean(receptor[r][c] for r in RECEPTORS) for c in bounds
        }
        results[arm] = {}
        for state, values in receptor.items():
            labels = [bounds[c][0] for c in point_ids]
            preds = [values[c] for c in point_ids]
            correct = {
                (a, b): 1
                if values[a] < values[b]
                else 0
                if values[a] > values[b]
                else 0.5
                for a, b in edges
            }
            flags[arm, state] = correct
            results[arm][state] = {
                "score_medians_or_mean": {str(c): float(v) for c, v in values.items()},
                "correct": sum(correct.values()),
                "pairs": len(edges),
                "spearman": statistics.correlation(ranks(labels), ranks(preds)),
                "kendall_tau_b": kendall(labels, preds),
            }
    changes = {}
    for rec in (*RECEPTORS, "receptor_mean"):
        base = flags["parent", rec]
        fixed = flags["fixed", rec]
        changes[rec] = {
            "gained": [k for k in edges if fixed[k] > base[k]],
            "lost": [k for k in edges if fixed[k] < base[k]],
        }
    keys = list(itertools.product(bounds, RECEPTORS, SEEDS))
    return {
        "n_raw_new_pools": 96,
        "n_official_archived_pools": 48,
        "interval_valued_compounds": [10],
        "aggregate": results,
        "parent_to_fixed_changes": changes,
        "changed_top_coordinate_sets": sum(
            coords["parent", *k][0] != coords["fixed", *k][0] for k in keys
        ),
        "changed_ordered_coordinate_pools": sum(
            coords["parent", *k] != coords["fixed", *k] for k in keys
        ),
        "any_changed_candidates_ignoring_order": sum(
            set(coords["parent", *k]) != set(coords["fixed", *k]) for k in keys
        ),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "cubic_grid", "out"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    result = check(a.work, a.cubic_grid)
    a.out.write_text(json.dumps(result, indent=2) + "\n")
