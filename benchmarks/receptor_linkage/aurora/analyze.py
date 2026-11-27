"""Reconstruct state masses, check them by kinetics, and audit the released docking scores."""

import argparse
import csv
import hashlib
import itertools
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R = 0.00198720425864083


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, data):
    with path.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=data[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def equilibrium(edges):
    # At free [L] = 1 µM the numerical weights are relative to one free 'in' state.
    masses = {"in": 1.0}
    pending = edges[:]
    while pending:
        old = len(pending)
        for e in pending[:]:
            a, b = e["start"], e["finish"]
            ratio = float(e["forward"]) / float(e["reverse"])
            if a in masses:
                value = masses[a] * ratio
                if b in masses and not math.isclose(masses[b], value, rel_tol=1e-12):
                    raise ValueError("Cycle violates detailed balance")
                masses[b] = value
            elif b in masses:
                masses[a] = masses[b] / ratio
            else:
                continue
            pending.remove(e)
        if len(pending) == old:
            raise ValueError("Disconnected state graph")
    free = masses["in"] + masses["out"]
    bound = sum(v for s, v in masses.items() if ":" in s)
    return 1000 * free / bound, masses["out"] / free


def stationary_check(edges, kd_nm):
    # Independent check using time evolution at three ligand activities.
    states = sorted({e[k] for e in edges for k in ("start", "finish")})
    errors = []
    for factor in (0.1, 1.0, 10.0):
        concentration_um = kd_nm * factor / 1000
        generator = np.zeros((len(states), len(states)))
        for e in edges:
            a, b = [states.index(e[k]) for k in ("start", "finish")]
            generator[a, b] = float(e["forward"]) * (
                concentration_um if e["step"] == "binding" else 1
            )
            generator[b, a] = float(e["reverse"])
        generator[np.diag_indices(len(states))] = -generator.sum(axis=1)
        system = generator.T.copy()
        system[-1] = 1
        rhs = np.zeros(len(states))
        rhs[-1] = 1
        probabilities = np.linalg.solve(system, rhs)
        occupied = sum(p for p, state in zip(probabilities, states) if ":" in state)
        errors.append(abs(occupied - factor / (1 + factor)))
    if max(errors) > 1e-9:
        raise ArithmeticError(errors)
    return max(errors)


def kinetics(out):
    result = []
    data = rows(HERE / "rates.csv")
    for name, subset in itertools.groupby(data, key=lambda e: e["compound"]):
        edges = list(subset)
        variants = [("printed_scheme", edges)]
        if name == "at9283":
            # Both 25 C ligands see the same apo direction as the independent Fig. 5 supplement 2.
            variants.append(
                (
                    "same_apo_as_5I",
                    [
                        dict(e, forward=e["reverse"], reverse=e["forward"])
                        if e["step"] == "apo"
                        else e
                        for e in edges
                    ],
                )
            )
        for variant, network in variants:
            kd, p_out = equilibrium(network)
            encounter, _ = equilibrium(
                [e for e in network if e["step"] != "relaxation"]
            )
            err = stationary_check(network, kd)
            experiment = float(edges[0]["wave_nM"])
            temperature = 273.15 + float(edges[0]["kinetics_C"])
            result.append(
                {
                    "compound": name,
                    "variant": variant,
                    "apo_out_fraction": p_out,
                    "encounter_nM": encounter,
                    "full_nM": kd,
                    "wave_nM": experiment,
                    "wave_sd_nM": edges[0]["wave_sd_nM"],
                    "relaxation_factor": encounter / kd,
                    "relaxation_dG_kcal_mol": R
                    * temperature
                    * math.log(kd / encounter),
                    "stationary_max_abs_error": err,
                }
            )
    write(out / "kinetic_reconstruction.csv", result)


def audit_docking(work, out):
    ledger = rows(work / "runs.csv")
    groups = defaultdict(list)
    medians = {}
    costs = []
    expected = set(
        itertools.product(
            ("official", "pro", "double"),
            ("2J4Z", "2X6E"),
            ("imatinib", "danusertib", "at9283"),
            (60913, 60943, 61211, 61231),
        )
    )
    for row in ledger:
        method, receptor, compound, seed = (
            row[k] for k in ("method", "receptor", "compound", "seed")
        )
        expected.remove((method, receptor, compound, int(seed)))
        stem = f"{method}_{receptor}_{compound}_{seed}"
        data = (work / f"pools/{stem}.pdbqt").read_bytes()
        scores = [
            float(line.split()[3])
            for line in data.splitlines()
            if line.startswith(b"REMARK VINA RESULT:")
        ]
        if (
            hashlib.sha256(data).hexdigest() != row["pool_sha256"]
            or len(scores) != int(row["modes"])
            or scores[0] != float(row["score"])
        ):
            raise ValueError(stem)
        log = (work / f"pools/{stem}.log").read_text()
        evals = [
            int(x)
            for x in re.findall(
                r"(?:hunt_evals|second_evals|explicit_evals)=(\d+)", log
            )
        ]
        row["reported_evals"] = sum(evals) if evals else ""
        groups[method, receptor, compound].append(scores[0])
    if expected:
        raise ValueError(expected)
    for key, scores in sorted(groups.items()):
        medians[key] = statistics.median(scores)
    stats, differences = [], []
    for (m, r, name), val in medians.items():
        stats.append(
            {
                "method": m,
                "receptor": r,
                "compound": name,
                "median_top_score": val,
                "minimum": min(groups[m, r, name]),
                "maximum": max(groups[m, r, name]),
            }
        )
    measured = {r["compound"]: float(r["wave_nM"]) for r in rows(HERE / "rates.csv")}
    for m, r in itertools.product(("official", "pro", "double"), ("2J4Z", "2X6E")):
        for a, b in itertools.combinations(("danusertib", "at9283", "imatinib"), 2):
            differences.append(
                {
                    "method": m,
                    "receptor": r,
                    "ligand_A": a,
                    "ligand_B": b,
                    "median_score_B_minus_A": medians[m, r, b] - medians[m, r, a],
                    "observed_kd_B_over_A": measured[b] / measured[a],
                    "same_running_buffer": (b != "imatinib"),
                }
            )
    for method in ("official", "pro", "double"):
        cases = [r for r in ledger if r["method"] == method]
        costs.append(
            {
                "method": method,
                "runs": len(cases),
                "modes": sum(int(r["modes"]) for r in cases),
                "wall_s": sum(float(r["wall_s"]) for r in cases),
                "cpu_s": sum(float(r["cpu_s"]) for r in cases),
                "reported_hunt_second_explicit_evals": sum(
                    r["reported_evals"] for r in cases
                )
                if all(r["reported_evals"] != "" for r in cases)
                else "",
            }
        )
    write(out / "scores.csv", stats)
    write(out / "differences.csv", differences)
    write(out / "costs.csv", costs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    a.output.mkdir(exist_ok=True, parents=True)
    kinetics(a.output)
    audit_docking(a.work, a.output)
