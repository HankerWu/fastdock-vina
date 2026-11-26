#!/usr/bin/env python3
"""Compare two ligands across matched rigid states with unknown apo populations.

The input energies must use a common unbound reference for each ligand. The
bounds quantify uncertainty in state weights conditional on these input
scores, not errors of the scoring function or unsampled states.
"""

import argparse
import csv
import itertools
import math
import sys
from collections import defaultdict
from statistics import median

R = 0.00198720425864083


def pair_sensitivity(first, second, kelvin=298.15):
    """Bound G(second)-G(first); positive values favor the first ligand."""
    if len(first) != len(second) or len(first) == 0:
        raise ValueError("Both ligands need the same nonempty set of states")
    if not all(math.isfinite(v) for v in [*first, *second, kelvin]) or kelvin <= 0:
        raise ValueError("Expected finite energies and positive temperature")
    gaps = [b - a for a, b in zip(first, second)]
    rt = R * kelvin

    def uniform_energy(values):
        center = min(values)
        return center - rt * math.log(
            sum(math.exp(-(v - center) / rt) for v in values) / len(values)
        )

    result = {
        "gap_min": min(gaps),
        "gap_max": max(gaps),
        "gap_equal_population": uniform_energy(second) - uniform_energy(first),
        "direction": "first"
        if min(gaps) > 0
        else "second"
        if max(gaps) < 0
        else "ambiguous",
        "p_first_state_crossover": None,
    }
    if len(first) == 2 and gaps[0] * gaps[1] < 0:
        # D_i = exp(-A_i/RT) - exp(-B_i/RT). Evaluating each magnitude
        # in log space avoids cancellation for small score gaps.
        def log_abs_difference(a, b):
            return -min(a, b) / rt + math.log(-math.expm1(-abs(b - a) / rt))

        h0, h1 = (log_abs_difference(a, b) for a, b in zip(first, second))
        if h0 >= h1:
            ratio = math.exp(h1 - h0)
            p = ratio / (1 + ratio)
        else:
            p = 1 / (1 + math.exp(h0 - h1))
        result["p_first_state_crossover"] = p
    return result


def run(source, dest, kelvin):
    rows = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for row in csv.DictReader(source):
        target, ligand, state, replicate = (
            row[key] for key in ("target", "ligand", "state", "replicate")
        )
        if not all((target, ligand, state, replicate)):
            raise ValueError("Empty target, ligand, state or replicate")
        values = rows[target][ligand][state]
        if replicate in values:
            raise ValueError(
                f"Duplicate observation: {target}, {ligand}, {state}, {replicate}"
            )
        value = float(row["energy"])
        if not math.isfinite(value):
            raise ValueError("Nonfinite energy")
        values[replicate] = value
    writer = csv.DictWriter(
        dest,
        fieldnames=[
            "target",
            "first",
            "second",
            "states",
            "gap_min",
            "gap_max",
            "gap_equal_population",
            "direction",
            "p_first_state_crossover",
        ],
    )
    writer.writeheader()
    for target, ligands in sorted(rows.items()):
        for first, second in itertools.combinations(sorted(ligands), 2):
            a, b = ligands[first], ligands[second]
            if a.keys() != b.keys() or any(a[s].keys() != b[s].keys() for s in a):
                raise ValueError(
                    f"Unmatched state or replicate: {target}, {first}, {second}"
                )
            states = sorted(a)
            # Seeds repeat a search; they receive no statistical-mechanical weight.
            features = pair_sensitivity(
                [median(a[s].values()) for s in states],
                [median(b[s].values()) for s in states],
                kelvin,
            )
            writer.writerow(
                dict(
                    target=target,
                    first=first,
                    second=second,
                    states="|".join(states),
                    **features,
                )
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "observations",
        type=argparse.FileType("r"),
        help="CSV columns target,ligand,state,replicate,energy",
    )
    parser.add_argument("--kelvin", type=float, default=298.15)
    args = parser.parse_args()
    run(args.observations, sys.stdout, args.kelvin)
