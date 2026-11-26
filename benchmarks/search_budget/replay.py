#!/usr/bin/env python3
"""Replay archived native docking pools and same-assay Kd score aggregations."""
import argparse
import csv
import itertools
import json
from collections import defaultdict
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
POSE_FILES = (
    "benchmarks/cubic_grid/human_panel/modes_native_24000.csv",
    "benchmarks/torsion_proposals/results/modes-native-12000-bridge.csv",
)
AFFINITY = "benchmarks/cubic_grid/affinity"


def rows(root, relative):
    with (root / relative).open(newline="") as f:
        yield from csv.DictReader(f)


def pose_replay(root, relative):
    grouped = defaultdict(lambda: defaultdict(list))
    for row in rows(root, relative):
        if row["kind"].lower() != "cross":
            continue
        key = (row["target"], row.get("query") or row["ligand"], row["receiver"])
        grouped[key][int(row["seed"])].append(
            (int(row.get("rank") or row["model"]), Decimal(row["score"]), float(row["rmsd"]))
        )
    panel = {}
    for key, seeds in grouped.items():
        if len(seeds) != 3:
            raise ValueError(f"{key}: expected three seeds")
        records = []
        for seed, pool in sorted(seeds.items()):
            ordered = sorted(pool)
            if ordered[0][0] != 1 or len({m[0] for m in pool}) != len(pool):
                raise ValueError(f"{key}/{seed}: invalid ranks")
            records.append((ordered[0][1], ordered[0][2] <= 2, any(m[2] <= 2 for m in pool), seed))
        panel[key] = records
    n = len(panel)
    if n == 0:
        raise ValueError(f"{relative}: no cross-docking data")

    def outcome(p):
        best = min(p, key=lambda r: (r[0], r[3]))
        return int(best[1]), int(any(r[2] for r in p))

    single = [sum(outcome([r])[i] for p in panel.values() for r in p) for i in range(2)]
    all_three = [sum(outcome(p)[i] for p in panel.values()) for i in range(2)]
    two = [0, 0]
    adaptive = [0, 0]
    uniform = [0.0, 0.0]
    k = n // 2
    for held_out in range(3):
        base = {key: [r for i, r in enumerate(p) if i != held_out] for key, p in panel.items()}
        chosen = set(sorted(base, key=lambda key: (-abs(base[key][0][0] - base[key][1][0]), key))[:k])
        for key, pair in base.items():
            b, a = outcome(pair), outcome(panel[key])
            for i in range(2):
                two[i] += b[i]
                adaptive[i] += a[i] if key in chosen else b[i]
                uniform[i] += b[i] + (k / n) * (a[i] - b[i])
    return dict(file=relative, unique_directions=n, unique_targets=len({key[0] for key in panel}),
                seed_runs=3*n, one=single, two=two, three=all_three, extra_runs=3*k,
                adaptive=adaptive, random_expected=uniform)


def affinity_replay(root):
    panel = {}
    for row in rows(root, AFFINITY + "/panel.csv"):
        key = (row["target"], row["cid"])
        if key in panel:
            raise ValueError(f"Duplicate compound {key}")
        panel[key] = (Decimal(row["kd_nm"]), int(row["heavy"]))
    scores = defaultdict(dict)
    for row in rows(root, AFFINITY + "/cells.csv"):
        key = (row["target"], row["cid"], row["arm"])
        seed = int(row["seed"])
        if seed in scores[key]:
            raise ValueError(f"Duplicate score {key}/{seed}")
        scores[key][seed] = Decimal(row["score"])
    result = {}
    for target in sorted({k[0] for k in panel}):
        ids = sorted(k[1] for k in panel if k[0] == target)
        pairs = []
        for a, b in itertools.combinations(ids, 2):
            if panel[(target,a)][0] > panel[(target,b)][0]:
                a,b = b,a
            if panel[(target,b)][0] >= 3 * panel[(target,a)][0]:
                pairs.append((a,b))

        def credit(a,b):
            return 1 if a < b else 0 if a > b else 0.5

        arms = {}
        for method in ("native", "cubic"):
            data = {cid: list(scores[(target,cid,method)].values()) for cid in ids}
            if any(len(x) != 2 for x in data.values()):
                raise ValueError(f"{target}/{method}: expected two scores per compound")
            arms[method] = {
                name: sum(credit(fn(data[a]),fn(data[b])) for a,b in pairs)
                for name,fn in (("min",min),("mean",lambda x: sum(x)/len(x)),("max",max))
            }
        arms["size_only"] = sum(credit(-panel[(target,a)][1],-panel[(target,b)][1]) for a,b in pairs)
        result[target] = dict(pairs=len(pairs), compounds=len(ids), **arms)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    print(json.dumps({"pose": [pose_replay(args.root, path) for path in POSE_FILES],
                      "affinity": affinity_replay(args.root)}, indent=2))


if __name__ == "__main__":
    main()
