"""Independent verification from author XLS, raw PDBQT and a direct p grid.

Counts jointly reproducible signs after collapsing repeated search seeds,
checks active compounds by an independent published data release and
computes paired bootstrap sensitivity to resampling whole ligands.
"""

import argparse
import csv
from collections import defaultdict
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import xlrd


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def audit(work, roots):
    panel = read(work / "panel.csv")
    measurements = xlrd.open_workbook(work / "raw/table4.xls").sheet_by_index(0)
    annotation = xlrd.open_workbook(work / "raw/table1.xls").sheet_by_index(0)
    names = list(json.loads((work / "raw/ligands_iso.txt").read_text()))
    proteins = list(json.loads((work / "raw/proteins.txt").read_text()))
    matrix = np.loadtxt(
        work / "raw/drug-target_interaction_affinities_Kd__Davis_et_al.2011v1.txt"
    )
    kd = {}
    for r in panel:
        target, drug = r["target"], r["cid"]
        row = measurements.col_values(2).index(target)
        col = measurements.row_values(0).index(r["compound"])
        value = float(measurements.cell_value(row, col))
        assert measurements.cell_value(row, 0) == r["accession"]
        assert (
            value
            == float(r["kd_nm"])
            == matrix[names.index(drug), proteins.index(target)]
        )
        assert annotation.cell_value(annotation.col_values(2).index(target), 3) == "NO"
        kd[target, drug] = value

    states = {
        "CDK2": ["1HCK", "1G5S"],
        "GSK3B": ["1Q3D", "1Q5K"],
        "SRC": ["1Y57", "2SRC"],
    }
    per_target = {}
    for target in states:
        compounds = [c for t, c in kd if t == target]
        per_target[target] = [
            (a, b)
            for a, b in itertools.permutations(compounds, 2)
            if kd[target, b] >= 3 * kd[target, a]
        ]
    destinations = {
        "native": ("pools", "results"),
        "mixed": ("mixed-pools", "mixed-results"),
        "soft": ("soft-pools", "soft-results"),
    }
    scoring = ["native", "fixed"]
    summary = []
    edge_details = []
    comparison = []
    for cohort, cohort_roots in roots.items():
        scores = defaultdict(lambda: defaultdict(list))
        for root in cohort_roots:
            for arm, (folder, result) in destinations.items():
                rows = read(root / result / "cells.csv")
                for row in rows:
                    target, cid, state, seed = (
                        row[k] for k in ("target", "ligand", "state", "replicate")
                    )
                    pool = root / f"{folder}/{state}_{cid}_{seed}.pdbqt"
                    blob = pool.read_bytes()
                    assert hashlib.sha256(blob).hexdigest() == row["pool_sha256"]
                    text = blob.decode().split("ENDMDL")[0]
                    native = float(text.split("VINA RESULT:")[1].split()[0])
                    raw = float(text.split("INTER + INTRA:")[1].split()[0])
                    assert native == float(row["native"])
                    fixed = float(row["normalization"]) * (
                        raw - float(row["prepared_unbound"])
                    )
                    assert abs(fixed - float(row["energy"])) < 1.0e-9
                    for scoring_name, value in zip(scoring, [native, fixed]):
                        scores[arm, target, scoring_name, cid][state].append(
                            (seed, value)
                        )
        for target, edges in per_target.items():
            mols = sorted({c for edge in edges for c in edge})
            for score in scoring:
                means = {}
                for arm in destinations:
                    for i, cid in enumerate(mols):
                        reps = [
                            scores[arm, target, score, cid][s] for s in states[target]
                        ]
                        assert set(r[0] for r in reps[0]) == set(r[0] for r in reps[1])
                        means[arm, cid] = np.array(
                            [np.median([r[1] for r in v]) for v in reps]
                        )
                    stats = defaultdict(list)
                    for a, b in edges:
                        ea, eb = (means[arm, c] for c in (a, b))
                        # Direct association proxy integration in ordinary space.
                        ka, kb = (
                            np.exp(-e / (298.15 * 0.00198720425864083))
                            for e in (ea, eb)
                        )
                        direct = (
                            math.log(ka.sum() / kb.sum()) * 298.15 * 0.00198720425864083
                        )
                        robust = min(eb - ea) > 0 or max(eb - ea) < 0
                        if len(cohort_roots) == 1:
                            published = next(
                                row
                                for row in read(
                                    cohort_roots[0] / destinations[arm][1] / "pairs.csv"
                                )
                                if row["target"] == target
                                and row["strong_cid"] == a
                                and row["weak_cid"] == b
                            )
                            assert (
                                abs(
                                    float(published[score + "_gap_equal_population"])
                                    - direct
                                )
                                < 1.0e-10
                            )
                            if published[score + "_p_first_state_crossover"]:
                                w = float(published[score + "_p_first_state_crossover"])
                                assert (
                                    abs(
                                        np.log(
                                            np.dot([w, 1 - w], ka)
                                            / np.dot([w, 1 - w], kb)
                                        )
                                    )
                                    < 1.0e-9
                                )
                        observation = {
                            "state0": float(eb[0] - ea[0]),
                            "state1": float(eb[1] - ea[1]),
                            "equal": direct,
                            "mean": float(np.mean(eb - ea)),
                        }
                        for name, delta in observation.items():
                            stats[name].append(
                                float(delta > 0) + 0.5 * float(delta == 0)
                            )
                        stats["robust"].append(robust)
                        stats["robust_correct"].append(robust and direct > 0)
                        edge_details.append(
                            dict(
                                cohort=cohort,
                                target=target,
                                score=score,
                                arm=arm,
                                strong=a,
                                weak=b,
                                kd_ratio=kd[target, b] / kd[target, a],
                                **observation,
                                robust=robust,
                            )
                        )
                    summary.append(
                        dict(
                            cohort=cohort,
                            target=target,
                            arm=arm,
                            score=score,
                            pairs=len(edges),
                            **{name: sum(v) for name, v in stats.items()},
                        )
                    )
                # Same ligand draw for native/mixed/soft: edges never independent.
                rng = np.random.default_rng(20261126)
                idx = {c: i for i, c in enumerate(mols)}
                draws = rng.multinomial(
                    len(mols), np.full(len(mols), 1 / len(mols)), size=10000
                )
                weight = np.stack(
                    [draws[:, idx[a]] * draws[:, idx[b]] for a, b in edges], axis=1
                )
                weight = weight[weight.sum(1) > 0]
                subset = [
                    r
                    for r in edge_details
                    if r["target"] == target
                    and r["cohort"] == cohort
                    and r["score"] == score
                ]
                for arm in ("mixed", "soft"):
                    for metric in ("equal", "mean"):
                        base = np.array(
                            [x[metric] for x in subset if x["arm"] == "native"]
                        )
                        candidate = np.array(
                            [x[metric] for x in subset if x["arm"] == arm]
                        )
                        signs = lambda e: np.sign(e) / 2 + 0.5
                        delta = signs(candidate) - signs(base)
                        boot = (weight @ delta) / weight.sum(axis=1)
                        lo, hi = np.quantile(boot, [0.025, 0.975])
                        comparison.append(
                            dict(
                                cohort=cohort,
                                target=target,
                                arm=arm,
                                score=score,
                                metric=metric,
                                pairs=len(edges),
                                gained=sum(delta > 0),
                                lost=sum(delta < 0),
                                net=sum(delta),
                                delta_pp=100 * delta.mean(),
                                lower_pp=100 * lo,
                                upper_pp=100 * hi,
                                draws=len(boot),
                            )
                        )
    write(work / "results/independent-summaries.csv", summary)
    write(work / "results/independent-contrasts.csv", comparison)
    write(work / "results/independent-edges.csv", edge_details)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--replicate", action="store_true")
    args = p.parse_args()
    roots = {"primary": [args.work]}
    if args.replicate:
        roots.update(
            {
                "replication": [args.work / "independent-seeds"],
                "all": [args.work, args.work / "independent-seeds"],
            }
        )
    audit(args.work, roots)
