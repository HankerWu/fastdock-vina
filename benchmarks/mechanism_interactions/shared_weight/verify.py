"""Independent direct-energy checks on all crossings and the open intervals."""
import argparse
import csv
import json
import sys
from bisect import bisect_left
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from statistics import fmean

HERE = Path(__file__).resolve().parent


def partition(pool, slack):
    cuts = {Fraction(0)}
    for a, b in combinations(pool, 2):
        if a.term == b.term:
            continue
        for epsilon in {slack, -slack}:
            crossing = (a.dry - b.dry + epsilon) / (b.term - a.term)
            if crossing >= 0:
                cuts.add(crossing)
    cuts = sorted(cuts)
    def loss(w):
        scores = [m.dry + w * m.term for m in pool]
        cutoff = min(scores) + slack
        return min(m.rmsd for m, score in zip(pool, scores) if score <= cutoff)

    # Evaluate every exact crossing, every open cell, and the unbounded tail.
    at = [loss(x) for x in cuts]
    right = [loss((a+b)/2) for a, b in zip(cuts, cuts[1:])] + [loss(cuts[-1]+1)]
    active = [i for i in range(len(cuts)) if i == 0 or
              len({right[i-1], at[i], right[i]}) != 1]
    return [[seq[i] for i in active] for seq in (cuts, at, right)]


def add_curves(parts, grid):
    cuts, at, right = parts
    res = []
    for w in grid:
        i = bisect_left(cuts, w)
        res.append(at[i] if i < len(cuts) and cuts[i] == w else right[i-1])
    return res


def check(pools, summaries, slack):
    parts = {k: partition(p, slack) for k, p in pools.items()}
    cuts = sorted(set(x for c, _, _ in parts.values() for x in c))
    weights = cuts + [(a+b)/2 for a, b in zip(cuts, cuts[1:])] + [cuts[-1]+1]
    curves = {k: add_curves(v, weights) for k, v in parts.items()}
    for row in summaries:
        group = row["group"].split("+")
        keys = [k for k in pools if k[0] == row["experiment"] and
                (row["panel"] == "all_historical" and k[1] != "gla" or k[1] == row["panel"])
                and (group == ["all_targets"] or k[2] in group)
                and ("search" not in row or k[4] == row["search"])]
        targets = sorted({k[2] for k in keys})
        assert len(keys) == int(row["pools"]) and sum(len(pools[k]) for k in keys) == int(row["modes"])
        per_target = []
        solo = []
        oracle = []
        dry = []
        unit = []
        for t in targets:
            subset = [k for k in keys if k[2] == t]
            per_target.append([fmean(x) for x in zip(*(curves[k] for k in subset))])
            solo.append(fmean(min(curves[k]) for k in subset))
            oracle.append(fmean(min(m.rmsd for m in pools[k]) for k in subset))
            for lam, out in ((0, dry), (1, unit)):
                out.append(fmean(min(pools[k], key=lambda m: (m.dry + lam*m.term, m.number)).rmsd for k in subset))
        for col, value in [
            ("pool_specific_bound_A", fmean(solo)), ("geometric_oracle_A", fmean(oracle)),
            ("target_specific_bound_A", fmean(min(v) for v in per_target)),
            ("shared_bound_A", min(fmean(x) for x in zip(*per_target))),
            ("dry_actual_A", fmean(dry)), ("unit_actual_A", fmean(unit)),
        ]:
            assert abs(float(row[col]) - value) < 1e-11, (row, col, value)
    return len(weights)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE.parent / "attainability")
    parser.add_argument("--results", type=Path, default=HERE / "results")
    args = parser.parse_args()
    sys.path.insert(0, str(args.source.resolve()))
    from analyze import PRINT_SLACK, read_pools
    pools = read_pools() | read_pools(gla=True)
    rows = []
    for name in ["summary.csv", "target_pairs.csv", "search_arms.csv"]:
        with (args.results / name).open() as f:
            rows.extend(csv.DictReader(f))
    checks = {}
    for label, slack in [("nominal", Fraction(0)), ("print_interval", PRINT_SLACK)]:
        selected = [row for row in rows if row["score_precision"] == label]
        checks[label] = dict(checked_rows=len(selected), states_checked=check(pools, selected, slack))
    print(json.dumps(dict(pools=len(pools), modes=sum(map(len, pools.values())), checks=checks), indent=2))


if __name__ == "__main__":
    main()
