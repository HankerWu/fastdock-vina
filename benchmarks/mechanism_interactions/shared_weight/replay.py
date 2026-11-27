"""Lower bounds on published-pool RMSD when pools share one reranking coefficient."""
import argparse
import csv
import hashlib
import itertools
import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from statistics import mean

HERE = Path(__file__).resolve().parent


def envelope(pool, slack):
    """Intervals where score(i) is within slack of the exact pool minimum."""
    intervals = []
    for candidate in pool:
        lower, upper = Fraction(0), None
        for rival in pool:
            slope = candidate.term - rival.term
            gap = rival.dry - candidate.dry + slack
            if not slope:
                if gap < 0:
                    break
            elif slope < 0:
                lower = max(lower, gap / slope)
            else:
                cut = gap / slope
                upper = cut if upper is None else min(upper, cut)
            if upper is not None and upper < lower:
                break
        else:
            intervals.append((candidate, lower, upper))
    return intervals


def loss_curves(pools, slack):
    intervals = {key: envelope(pool, slack) for key, pool in pools.items()}
    # With label-favourable ties the minimum over any open interval is also
    # available at its closed left endpoint, simultaneously in every pool.
    # Thus computing all finite starts is sufficient, including for slack>0.
    weights = sorted({lo for partition in intervals.values() for _, lo, _ in partition})
    return weights, {key: [min(m.rmsd for m, lo, hi in parts
                               if lo <= w and (hi is None or w <= hi)) for w in weights]
                     for key, parts in intervals.items()}


def table_row(experiment, panel, group, pools, curves, target_keys):
    targets = sorted(target_keys)
    keys = [k for target in targets for k in target_keys[target]]
    # Equal targets, equal published pools within target. Curves at the same
    # index refer to the same coefficient for all targets.
    target_curves = [[mean(v) for v in zip(*(curves[k] for k in target_keys[t]))]
                     for t in targets]
    def average(f):
        return mean(mean(f(k) for k in target_keys[t]) for t in targets)
    oracle = average(lambda k: min(m.rmsd for m in pools[k]))
    pool_specific = average(lambda k: min(curves[k]))
    target_specific = mean(min(row) for row in target_curves)
    common = min(mean(v) for v in zip(*target_curves))
    dry = average(lambda k: min(pools[k], key=lambda m: (m.dry, m.number)).rmsd)
    unit = average(lambda k: min(pools[k], key=lambda m: (m.dry + m.term, m.number)).rmsd)
    bounds = [oracle, pool_specific, target_specific, common, min(dry, unit)]
    assert all(a <= b + 1e-12 for a, b in zip(bounds, bounds[1:])), bounds
    return dict(experiment=experiment, panel=panel, group=group, targets=len(targets),
                pools=len(keys), modes=sum(len(pools[k]) for k in keys),
                geometric_oracle_A=oracle, pool_specific_bound_A=pool_specific,
                target_specific_bound_A=target_specific, shared_bound_A=common,
                dry_actual_A=dry, unit_actual_A=unit)


def rows(pools, curves):
    by_target = defaultdict(list)
    for key in pools:
        by_target[key[2]].append(key)
    keys = list(pools)
    experiment, panel = keys[0][:2]
    if len({k[1] for k in keys}) > 1:
        panel = "all_historical"
    summary = [table_row(experiment, panel, "all_targets", pools, curves, by_target)]
    pairs = []
    for t in sorted(by_target):
        summary.append(table_row(experiment, panel, t, pools, curves, {t: by_target[t]}))
    for a, b in itertools.combinations(sorted(by_target), 2):
        pairs.append(table_row(experiment, panel, a + "+" + b, pools, curves,
                               {a: by_target[a], b: by_target[b]}))
    return summary, pairs


def panels(read_pools):
    historical = read_pools()
    for experiment in sorted({k[0] for k in historical}):
        subset = {k: v for k, v in historical.items() if k[0] == experiment}
        yield subset
        for panel in sorted({k[1] for k in subset}):
            yield {k: v for k, v in subset.items() if k[1] == panel}
    yield read_pools(gla=True)


def write(path, data):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=data[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE.parent / "attainability")
    parser.add_argument("--output", type=Path, default=HERE / "results")
    args = parser.parse_args()
    sys.path.insert(0, str(args.source.resolve()))
    from analyze import PRINT_SLACK, read_pools, sources

    summaries, pairs, search_arms = [], [], []
    for pools in panels(read_pools):
        for name, slack in [("nominal", Fraction(0)), ("print_interval", PRINT_SLACK)]:
            _, curves = loss_curves(pools, slack)
            s, p = rows(pools, curves)
            for arm in sorted({k[4] for k in pools}):
                subset = {k: v for k, v in pools.items() if k[4] == arm}
                by_target = defaultdict(list)
                for key in subset:
                    by_target[key[2]].append(key)
                record = table_row(s[0]["experiment"], s[0]["panel"], "all_targets",
                                   subset, curves, by_target)
                search_arms.append(dict(score_precision=name, search=arm, **record))
            for collection, batch in [(summaries, s), (pairs, p)]:
                collection.extend(dict(score_precision=name, **row) for row in batch)
    args.output.mkdir(parents=True, exist_ok=True)
    write(args.output / "summary.csv", summaries)
    write(args.output / "target_pairs.csv", pairs)
    write(args.output / "search_arms.csv", search_arms)
    paths = [p for gla in [False, True] for _, _, p in sources(gla)]
    hashes = {str(p.relative_to(args.source.resolve().parents[2])):
              hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (args.output / "inputs.json").write_text(json.dumps(hashes, indent=2) + "\n")


if __name__ == "__main__":
    main()
