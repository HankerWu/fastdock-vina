"""Exact lower-envelope diagnostics on the two archived factorial pose panels."""

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from statistics import mean, median

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
KEY = ("experiment", "panel", "target", "cap", "search", "query", "receiver", "seed")
GROUPS = KEY[:5]
PRINT_SLACK = Fraction(1, 1000)  # Difference between opposite half-unit endpoints


@dataclass(frozen=True)
class Mode:
    number: int
    dry: Fraction
    term: Fraction
    rmsd: float


def sources(gla=False):
    if gla:
        yield "ch_pi", "gla", BENCH / "ch_pi_search/gla_transfer/results/modes.csv"
        return
    for panel in ("gaa", "historical"):
        yield "ch_pi", panel, BENCH / "ch_pi_search" / "results" / panel / "modes.csv"
    yield "water", None, BENCH / "soft_water_interaction/results/poses.csv"


def read_pools(gla=False):
    pools = defaultdict(list)
    for experiment, panel, path in sources(gla):
        with path.open() as f:
            for row in csv.DictReader(f):
                if experiment == "water":
                    target, subpanel = row["target"], row["panel"]
                    term, dry = Fraction(row["water_penalty"]), Fraction(row["dry"])
                    receiver, number = row["receiver"], int(row["model"])
                else:
                    target, subpanel = row["panel"], panel
                    # Published Gaussian prefactor and within-ligand rotor scaling
                    term = Fraction(3, 10) * Fraction(row["raw_chpi"]) * Fraction(row["factor"])
                    dry = Fraction(row["native"])
                    receiver, number = row["receptor"], int(row["rank"])
                key = (experiment, subpanel, target, int(row["cap"]), row["search"],
                       row["query"], receiver, int(row["seed"]))
                pools[key].append(Mode(number, dry, term, float(row["rmsd"])))
    for key, modes in pools.items():
        modes.sort(key=lambda m: m.number)
        if [m.number for m in modes] != list(range(1, len(modes) + 1)):
            raise ValueError(f"Missing or repeated mode in {key}")
    return pools


def attainable(candidate, pool, slack=Fraction(0)):
    """Intersect every half-line E_i(w) <= E_j(w) + slack with finite w >= 0."""
    lower, upper = Fraction(0), None
    for rival in pool:
        delta_term = candidate.term - rival.term
        intercept_gap = rival.dry - candidate.dry + slack
        if delta_term == 0:
            if intercept_gap < 0:
                return False
        elif delta_term > 0:
            bound = intercept_gap / delta_term
            upper = bound if upper is None else min(upper, bound)
        else:
            lower = max(lower, intercept_gap / delta_term)
        if upper is not None and lower > upper:
            return False
    return True


def dominated(candidate, pool):
    return any(
        m.dry <= candidate.dry and m.term <= candidate.term
        and (m.dry < candidate.dry or m.term < candidate.term) for m in pool
    )


def summarize_pool(key, pool):
    identity = dict(zip(KEY, key))
    nominal = [m for m in pool if attainable(m, pool)]
    possible = [m for m in pool if attainable(m, pool, PRINT_SLACK)]
    for weight in (0, 1):
        best = min(m.dry + weight * m.term for m in pool)
        assert all(m in nominal for m in pool if m.dry + weight * m.term == best)
    assert set(nominal) <= set(possible)
    dry = min(pool, key=lambda m: (m.dry, m.number))
    unit = min(pool, key=lambda m: (m.dry + m.term, m.number))
    oracle = min(m.rmsd for m in pool)
    supported = min(m.rmsd for m in nominal)
    permissive = min(m.rmsd for m in possible)
    assert oracle <= permissive <= supported <= min(dry.rmsd, unit.rmsd)
    flags = []
    for m in pool:
        flags.append(dict(
            **identity, mode=m.number, oracle=int(m.rmsd == oracle),
            dry=int(m.number == dry.number), unit=int(m.number == unit.number),
            rmsd=m.rmsd, pareto_dominated=int(dominated(m, pool)),
            robustly_dominated=int(any(other.dry + PRINT_SLACK < m.dry
                                       and other.term <= m.term for other in pool)),
            attainable=int(m in nominal), rounding_possible=int(m in possible),
        ))
    oracle_classes = {row["pareto_dominated"] for row in flags if row["oracle"]}
    cell = dict(
        **identity, modes=len(pool), attainable_modes=len(nominal),
        possible_modes=len(possible), dry_mode=dry.number, unit_mode=unit.number,
        dry_rmsd=dry.rmsd, unit_rmsd=unit.rmsd, pool_oracle=oracle,
        attainable_oracle=supported, rounding_oracle=permissive,
        nominal_gap=supported-oracle, rounding_gap=permissive-oracle,
        oracle_unattainable=int(supported > oracle),
        oracle_excluded_after_rounding=int(permissive > oracle),
        all_oracles_dominated=int(oracle_classes == {1}),
    )
    return cell, flags


def quantile(values, p):
    values = sorted(values)
    i, remainder = divmod((len(values) - 1) * p, 100)
    return values[i] + (values[min(i+1, len(values)-1)] - values[i]) * remainder / 100


def subset_summary(rows):
    gaps = [r["nominal_gap"] for r in rows]
    return dict(
        pools=len(rows), modes=sum(r["modes"] for r in rows),
        oracle_unattainable=sum(r["oracle_unattainable"] for r in rows),
        excluded_dominated=sum(r["oracle_unattainable"] and r["all_oracles_dominated"] for r in rows),
        excluded_nondominated=sum(r["oracle_unattainable"] and not r["all_oracles_dominated"] for r in rows),
        excluded_after_rounding=sum(r["oracle_excluded_after_rounding"] for r in rows),
        mean_pool_oracle=mean(r["pool_oracle"] for r in rows),
        mean_attainable_oracle=mean(r["attainable_oracle"] for r in rows),
        mean_dry_rmsd=mean(r["dry_rmsd"] for r in rows),
        mean_unit_rmsd=mean(r["unit_rmsd"] for r in rows),
        mean_gap=mean(gaps), median_gap=median(gaps), p90_gap=quantile(gaps, 90),
        mean_rounding_gap=mean(r["rounding_gap"] for r in rows),
    )


def write_rows(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(out, gla=False):
    out.mkdir(parents=True, exist_ok=True)
    cells, modes = [], []
    for key, pool in sorted(read_pools(gla).items()):
        cell, flags = summarize_pool(key, pool)
        cells.append(cell)
        modes.extend(flags)
    groups = defaultdict(list)
    for r in cells:
        groups[tuple(r[k] for k in GROUPS)].append(r)
    summaries = [dict(**dict(zip(GROUPS, k)), **subset_summary(v)) for k, v in groups.items()]
    distributions = []
    for k, rows in groups.items():
        for metric in ("dry_rmsd", "unit_rmsd", "pool_oracle", "attainable_oracle", "nominal_gap"):
            values = [row[metric] for row in rows]
            distributions.append(dict(**dict(zip(GROUPS, k)), metric=metric, pools=len(values),
                                      mean=mean(values), q25=quantile(values, 25),
                                      median=median(values), q75=quantile(values, 75),
                                      p90=quantile(values, 90), p95=quantile(values, 95)))
    write_rows(out / "distributions.csv", distributions)
    write_rows(out / "modes.csv", modes)
    write_rows(out / "pools.csv", cells)
    write_rows(out / "strata.csv", summaries)
    inputs = {str(p.relative_to(BENCH.parent)): hashlib.sha256(p.read_bytes()).hexdigest()
              for _, _, p in sources(gla)}
    (out / "input_sha256.json").write_text(json.dumps(inputs, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--gla", action="store_true", help="Use only the separate GLA transfer panel")
    args = parser.parse_args()
    run(args.output or HERE / ("gla_results" if args.gla else "results"), args.gla)
