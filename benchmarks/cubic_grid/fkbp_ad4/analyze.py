"""Recompute FKBP12 series orders and all AD4 mode scores directly from maps and printed coordinates."""

import argparse
import csv
import itertools
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import spearmanr

from run_official import COMPOUNDS, RECEPTORS, SEEDS, TYPES, sha, table


def read_map(path):
    with path.open() as f:
        headers = [next(f).split() for _ in range(6)]
        spacing = float(headers[3][1])
        n = np.array([int(x) for x in headers[4][1:]])
        center = np.array([float(x) for x in headers[5][1:]])
        values = np.loadtxt(f).reshape(tuple(n + 1)[::-1]).T
    lo = center - n * spacing / 2
    axes = [x + spacing * np.arange(m + 1) for x, m in zip(lo, n)]
    return RegularGridInterpolator(axes, values, method="linear", bounds_error=True)


def calculate_components(lines, fields, *, operator="separate", uncapped_terms=False):
    """Reconstruct an in-box, rigid-receptor AD4 endpoint from printed atoms.

    Official 1.2.7 caps the three positive map values separately. The Pro
    combined operator caps the assembled signed atom field. Return the legacy
    official breakdown by default. For either operator, uncapped_terms=True
    separates the raw additive maps from the nonadditive cap remainder.
    """
    if operator not in ("separate", "combined"):
        raise ValueError(f"Unknown AD4 map operator: {operator}")
    atoms = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    coords = np.array([[float(a[30 + k:38 + k]) for k in (0, 8, 16)] for a in atoms])
    charge = np.array([float(a[70:76]) for a in atoms])
    kinds = np.array([a[77:].strip() for a in atoms])
    dof = int(next(line for line in lines if line.startswith("TORSDOF")).split()[1])

    def val(kind, xyz):
        field = fields[kind]
        lo = [axis[0] for axis in field.grid]
        hi = [axis[-1] for axis in field.grid]
        # Printed coordinates at a decimal box face can differ by one floating-point ulp.
        inside = np.clip(xyz, lo, hi)
        if np.max(np.abs(inside - xyz)) > 1e-8:
            raise ValueError("A written atom lies outside the scoring box")
        return field(inside)

    def cap(energy):
        capped = energy.copy()
        positive = capped > 0
        capped[positive] *= 1000 / (1000 + energy[positive])
        return capped

    typed = np.empty(len(atoms))
    for kind in sorted(set(kinds)):
        mask = kinds == kind
        typed[mask] = val(kind, coords[mask])
    electric = val("e", coords)
    desolv = val("d", coords)
    if operator == "separate" and not uncapped_terms:
        # Keep the released official component definitions and summation order.
        return dict(
            type_maps=sum(cap(typed[kinds == kind]).sum() for kind in sorted(set(kinds))),
            electrostatic=float(charge @ cap(electric)),
            charge_desolvation=float(abs(charge) @ cap(desolv)),
            torsion=0.2983 * dof,
        )

    raw = typed + charge * electric + abs(charge) * desolv
    if operator == "combined":
        final = cap(raw)
    else:
        final = cap(typed) + charge * cap(electric) + abs(charge) * cap(desolv)
    return dict(
        type_maps=typed.sum(),
        electrostatic=float(charge @ electric),
        charge_desolvation=float(abs(charge) @ desolv),
        cap_adjustment=(final - raw).sum(),
        torsion=0.2983 * dof,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "labels", "out"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    a.out.mkdir(exist_ok=True, parents=True)
    labels = {
        int(l["compound"]): (float(l["ki_low_nm"]), float(l["ki_high_nm"]))
        for l in csv.DictReader(a.labels.open())
    }
    assert set(labels) == set(COMPOUNDS)
    official = list(csv.DictReader((a.work / "cells.csv").open()))
    fingerprint = {
        (r["receiver"], int(r["compound"]), int(r["seed"]), r["arm"]): r
        for r in official
    }
    tasks = list(itertools.product(RECEPTORS, COMPOUNDS, SEEDS, ("ad4", "vina")))
    assert len(official) == len(tasks) and set(fingerprint) == set(tasks)
    read = {
        r: {
            t: read_map(a.work / "maps" / r / f"receptor.{t}.map")
            for t in (*TYPES, "e", "d")
        }
        for r in RECEPTORS
    }
    pools = defaultdict(list)
    audit = []
    torsion_sensitivity = defaultdict(list)
    for r, n, seed, arm in tasks:
        stem = f"{n}_{r}_{seed}_{arm}"
        path = a.work / "pools" / (stem + ".pdbqt")
        log = (a.work / "pools" / (stem + ".log")).read_text()
        assert f"Scoring function : {arm}" in log and f"(random seed: {seed})" in log
        assert sha(path) == fingerprint[(r, n, seed, arm)]["sha256"]
        modes = [
            b.strip().splitlines()
            for b in path.read_text().split("ENDMDL")
            if b.strip()
        ]
        assert len(modes) == int(fingerprint[(r, n, seed, arm)]["modes"])
        logged_scores = {}
        for line in log.splitlines():
            fields = line.split()
            if len(fields) == 4 and fields[0].isdigit():
                logged_scores[int(fields[0])] = float(fields[1])
        assert set(logged_scores) == set(range(1, len(modes) + 1))
        rescores = []
        for rank, lines in enumerate(modes, 1):
            assert lines[0].split() == ["MODEL", str(rank)]
            reported = float(
                next(l for l in lines if l.startswith("REMARK VINA RESULT:")).split()[3]
            )
            # The CLI table uses four significant digits, the pose file three decimal places.
            assert abs(logged_scores[rank] - reported) < 0.0051
            if rank == 1:
                assert reported == float(fingerprint[(r, n, seed, arm)]["score"])
                pools[(arm, r, n)].append(reported)
            comp = calculate_components(lines, read[r])
            energy = sum(comp.values())
            rescores.append(energy)
            audit.append(
                dict(
                    receiver=r,
                    compound=n,
                    seed=seed,
                    search_arm=arm,
                    mode=rank,
                    reported=reported,
                    ad4_from_maps=energy,
                    **comp,
                )
            )
        # The TORSDOF penalty is constant for all poses of the same ligand.
        if arm == "ad4":
            torsion_sensitivity[r, n].append(min(rescores) - comp["torsion"])
        if arm == "vina":
            pools[("ad4_on_vina_pool", r, n)].append(min(rescores))
    audit_err = [
        abs(x["ad4_from_maps"] - x["reported"])
        for x in audit
        if x["search_arm"] == "ad4"
    ]
    if max(audit_err) > 0.03:
        raise ValueError(f"AD4 score reconstruction error {max(audit_err)}")
    table(a.out / "mode_audit.csv", audit)
    methods = ("vina", "ad4", "ad4_on_vina_pool")
    mols = []
    scores = {}
    for method, n in itertools.product(methods, COMPOUNDS):
        rd = [float(np.median(pools[(method, r, n)])) for r in RECEPTORS]
        assert all(len(pools[(method, r, n)]) == 3 for r in RECEPTORS)
        scores[method, n] = sum(rd) / 2
        mols.append(
            dict(
                method=method,
                compound=n,
                median_1FKG=rd[0],
                median_1FKH=rd[1],
                score=scores[method, n],
            )
        )
    table(a.out / "scores.csv", mols)
    pairs = []
    for i, j in itertools.permutations(COMPOUNDS, 2):
        if labels[j][0] < 3 * labels[i][1]:
            continue
        row = dict(stronger=i, weaker=j, minimum_fold=labels[j][0] / labels[i][1])
        for method in methods:
            delta = scores[method, j] - scores[method, i]
            row[method + "_delta"] = delta
            row[method + "_credit"] = 0.5 if abs(delta) < 1e-10 else int(delta > 0)
        pairs.append(row)
    table(a.out / "pairs.csv", pairs)
    seven = [n for n, (lo, hi) in labels.items() if lo == hi]
    summary = {
        "n_compounds": len(labels),
        "n_unique_pairs": len(pairs),
        "n_runs": len(tasks),
        "n_modes": len(audit),
        "n_ad4_modes_reconstructed": len(audit_err),
        "max_ad4_map_rescore_error": max(audit_err),
        "ad4_gains": sum(r["ad4_credit"] > r["vina_credit"] for r in pairs),
        "ad4_losses": sum(r["ad4_credit"] < r["vina_credit"] for r in pairs),
        "methods": {
            m: {
                "correct": sum(row[m + "_credit"] for row in pairs),
                "rho_seven_points": float(
                    spearmanr(
                        [labels[n][0] for n in seven], [scores[m, n] for n in seven]
                    ).statistic
                ),
                "epimer_8_vs_10": scores[m, 10] - scores[m, 8],
                "epimer_9_vs_11": scores[m, 11] - scores[m, 9],
            }
            for m in methods
        },
    }
    no_torsion = {
        n: float(np.mean([np.median(torsion_sensitivity[r, n]) for r in RECEPTORS]))
        for n in COMPOUNDS
    }
    summary["exploratory_remove_torsion_from_ad4_pool"] = {
        "pairs_correct": sum(
            no_torsion[x["stronger"]] < no_torsion[x["weaker"]] for x in pairs
        ),
        "scores": no_torsion,
    }
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
