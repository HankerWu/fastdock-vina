"""Check envelope membership at line intersections and compare archived winners."""

import argparse
import csv
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path

from analyze import HERE, KEY, PRINT_SLACK, Mode, attainable, read_pools, sources


def intersections(lines):
    """All nonnegative abscissas at which two affine lines are equal, and zero."""
    points = {Fraction(0)}
    for (di, si), (dj, sj) in combinations(lines, 2):
        if si != sj:
            x = (dj - di) / (si - sj)
            if x >= 0:
                points.add(x)
    return points


def exposed_modes(lines):
    # Every winner interval starts at zero or a finite intersection. Ties at an
    # isolated point and lines coincident on the whole domain both count.
    supported = set()
    for x in intersections(lines):
        values = [intercept + x * slope for intercept, slope in lines]
        lowest = min(values)
        supported.update(i for i, y in enumerate(values) if y == lowest)
    return supported


def possible_with_rounding(i, lines):
    # Hold the other dry scores at their upper print limits, the tested one
    # at its lower limit. Translating all scores leaves it 0.001 lower.
    di, si = lines[i]
    di -= PRINT_SLACK
    # A nonempty feasible interval starts at zero or at a crossing of this
    # candidate (not necessarily two other lines) with a competitor.
    points = {Fraction(0)}
    rival_lines = [line for j, line in enumerate(lines) if j != i]
    for dj, sj in rival_lines:
        if sj != si:
            x = (dj - di) / (si - sj)
            if x >= 0:
                points.add(x)
    return any(all(di + si*x <= dj + sj*x for dj, sj in rival_lines) for x in points)


def check_limits():
    def check(pairs, expected):
        pool = [Mode(i, Fraction(d), Fraction(t), 0.0) for i, (d, t) in enumerate(pairs)]
        assert {m.number for m in pool if attainable(m, pool)} == expected
        assert exposed_modes([(m.dry, m.term) for m in pool]) == expected
    check([(3, 4)], {0})
    check([(0, 1), (0, 1), (1, 1)], {0, 1})
    check([(1, 1), (0, 2), (2, 0)], {0, 1, 2})  # single-point middle winner
    check([(2, 1), (0, 2), (3, 0)], {1, 2})      # unsupported nondominated
    check([(0, 2), (0, 1), (2, 0)], {0, 1, 2})  # winner only at w=0
    check([(0, 0), ('1e20', '-1e-30')], {0, 1})
    check([(0, 0), (-1, -1), (0, 1)], {1})
    rounded = [Mode(0, Fraction(0), Fraction(0), 0.0),
               Mode(1, PRINT_SLACK, Fraction(0), 1.0)]
    assert not attainable(rounded[1], rounded)
    assert attainable(rounded[1], rounded, PRINT_SLACK)
    assert possible_with_rounding(1, [(m.dry, m.term) for m in rounded])


def key(row):
    return tuple(int(row[k]) if k in ("cap", "seed") else row[k] for k in KEY)


def archive_cells(gla=False):
    for experiment, panel, mode_path in sources(gla):
        path = mode_path.with_name("cells.csv")
        with path.open() as f:
            for row in csv.DictReader(f):
                if experiment == "water":
                    target, subpanel, receiver = row["target"], row["panel"], row["receiver"]
                    weight, mode = int(row["ranker"] == "water"), row["model"]
                    rmsd, oracle = row["rmsd"], row["oracle"]
                else:
                    target, subpanel, receiver = row["panel"], panel, row["receptor"]
                    weight, mode = int(row["rerank"] == "field"), row["winner"]
                    rmsd, oracle = row["top1"], row["pool_min"]
                k = (experiment, subpanel, target, int(row["cap"]), row["search"],
                     row["query"], receiver, int(row["seed"]))
                yield k, weight, int(mode), float(rmsd), float(oracle)


def main(results, gla=False):
    check_limits()
    pools = read_pools(gla)
    with (results / "modes.csv").open() as f:
        mode_rows = {(key(r), int(r["mode"])): r for r in csv.DictReader(f)}
    with (results / "pools.csv").open() as f:
        cell_rows = {key(r): r for r in csv.DictReader(f)}
    assert len(mode_rows) == sum(map(len, pools.values())) and cell_rows.keys() == pools.keys()
    checked = set()
    for k, weight, mode, rmsd, oracle in archive_cells(gla):
        row = cell_rows[k]
        prefix = "unit" if weight else "dry"
        assert int(row[prefix + "_mode"]) == mode
        assert abs(float(row[prefix + "_rmsd"]) - rmsd) < 1e-12
        assert abs(float(row["pool_oracle"]) - oracle) < 1e-12
        assert (k, weight) not in checked
        checked.add((k, weight))
    assert len(checked) == 2 * len(pools)
    ambiguous_dry, ambiguous_unit = [], []
    changed_membership = 0
    for k, pool in pools.items():
        lines = [(m.dry, m.term) for m in pool]
        nominal = exposed_modes(lines)
        for i, m in enumerate(pool):
            row = mode_rows[(k, m.number)]
            possible = possible_with_rounding(i, lines)
            assert int(row["attainable"]) == (i in nominal)
            assert int(row["rounding_possible"]) == possible
            changed_membership += (i in nominal) != possible
        for w, ambiguous in [(0, ambiguous_dry), (1, ambiguous_unit)]:
            values = [d + w * t for d, t in lines]
            eligible = [m.number for m, v in zip(pool, values) if v <= min(values) + PRINT_SLACK]
            if len(eligible) > 1:
                ambiguous.append(dict(**dict(zip(KEY, k)), possible_modes=eligible))
    return dict(poses=len(mode_rows), pools=len(pools), archived_winners=len(checked),
                membership_changes_allowing_print_rounding=changed_membership,
                dry_readout_ambiguous=ambiguous_dry, unit_readout_ambiguous=ambiguous_unit)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path)
    parser.add_argument("--gla", action="store_true", help="Verify only the separate GLA transfer panel")
    args = parser.parse_args()
    result = main(args.results or HERE / ("gla_results" if args.gla else "results"), args.gla)
    print(json.dumps(result, indent=2))
