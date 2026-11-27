"""Independently reconstruct pKd, scores, compound-deletion effects and actual costs."""

import argparse
import hashlib
import itertools
import json
import re
from decimal import Decimal
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau, spearmanr
from transfer import COUNT, HERE, METHODS, SEED, calibrate, read, write

RECS = ("1KJL", "3T1L")
CONTRASTS = (
    ("parent", "short_parent"),
    ("parent", "retention"),
    ("parent", "recurrence"),
    ("short_parent", "retention"),
    ("short_parent", "recurrence"),
    ("retention", "recurrence"),
    ("parent", "full_retention"),
    ("parent", "full_recurrence"),
    ("full_retention", "full_recurrence"),
    ("full_retention", "retention"),
    ("full_recurrence", "recurrence"),
)


def digest(s):
    return hashlib.sha256(s.encode()).hexdigest()


def obtain(work):
    calibration = read(work / "calibration.csv")
    base_setup = json.loads((work / "setup.json").read_text())
    extra_setup = json.loads((work / "full_cap_setup.json").read_text())
    expected = hashlib.sha256(
        (HERE.parent / "affinity/cost.csv").read_bytes()
    ).hexdigest()
    assert (
        base_setup["calibration_sha256"]
        == extra_setup["calibration_sha256"]
        == expected
    )
    for method in ("retention", "recurrence"):
        assert (
            base_setup["binaries"][method] == extra_setup["binaries"]["full_" + method]
        )
    caps = {r["receiver"]: int(r["cap"]) for r in calibration}
    assert caps == {r["receiver"]: r["cap"] for r in calibrate()}
    data = {}
    rows = read(work / "runs.csv") + read(work / "full_cap_runs.csv")
    for r in rows:
        m, rec, c = r["arm"], r["receptor"], r["compound"]
        assert int(r["seed"]) == SEED and int(r["cap"]) == (
            12000 if m == "parent" or m.startswith("full_") else caps[rec]
        )
        prefix = work / "raw" / r["key"]
        text = prefix.with_suffix(".pdbqt").read_text()
        values = re.findall(r"^REMARK VINA RESULT:\s+(\S+)", text, re.MULTILINE)
        assert len(values) == int(r["modes"]) and values[0] == r["first_score"]
        first = text.split("ENDMDL", 1)[0]
        coords = "\n".join(
            l for l in first.splitlines() if l.startswith(("ATOM", "HETATM"))
        )
        log = prefix.with_suffix(".log").read_text()
        for k in COUNT:
            actual = re.findall(r"\b" + k + r"=(\d+)", log)
            assert actual == [r[k]]
        assert int(r["joint"]) == sum(int(r[k]) for k in COUNT[:3])
        explicit = int(re.search(r"\bexplicit_inputs=(\d+)", log)[1])
        decision = re.search(r"Deferred hunt:.*", log)
        fields = dict(re.findall(r"(\w+)=(\S+)", decision[0])) if decision else {}
        key = m, rec, c
        assert key not in data
        data[key] = {
            **r,
            "pose_sha256": digest(text),
            "first_xyz_sha256": digest(coords),
            "explicit_inputs": explicit,
            "deferred_replaced": fields.get("replaced", 0),
            "deferred_length": fields.get("run_length", ""),
            "deferred_skipped": fields.get("skipped", ""),
            "deferred_checks": fields.get("checks", ""),
        }
    assert len(data) == 564
    for rec, c in itertools.product(RECS, {r["compound"] for r in data.values()}):
        for m in ("retention", "recurrence", "full_retention", "full_recurrence"):
            reference = "parent" if m.startswith("full_") else "short_parent"
            for k in (
                "hunt_calls",
                "hunt_evals",
                "second_calls",
                "second_evals",
                "explicit_inputs",
            ):
                assert data[reference, rec, c][k] == data[m, rec, c][k]
    return data


def analyze(work, out):
    manifest = sorted(
        read(HERE.parents[1] / "pair_gradient/galectin_external/manifest.csv"),
        key=lambda x: int(x["compound"]),
    )
    n = len(manifest)
    assert n == 47
    ids = [x["compound"] for x in manifest]
    k = [Decimal(x["kd_uM"]) for x in manifest]
    pk = 6 - np.log10([float(x) for x in k])
    desc = [(int(x["heavy"]), int(x["charge"])) for x in manifest]
    data = obtain(work)
    estimates = []
    costs = []
    metrics = []
    flips = []
    deletions = []
    for rec in (*RECS, "mean"):
        group = RECS if rec == "mean" else (rec,)
        scores = {
            m: np.array(
                [
                    sum(
                        int(Decimal(data[m, s, c]["first_score"]) * 1000) for s in group
                    )
                    for c in ids
                ]
            )
            for m in (*METHODS, "full_retention", "full_recurrence")
        }
        for m in (*METHODS, "full_retention", "full_recurrence"):
            for i, c in enumerate(ids):
                estimates.append(
                    {
                        "receiver": rec,
                        "method": m,
                        "compound": c,
                        "Kd_uM": k[i],
                        "pKd": round(pk[i], 9),
                        "heavy": desc[i][0],
                        "formal_charge": desc[i][1],
                        "score_kcal_mol": scores[m][i] / (1000 * len(group)),
                    }
                )
            runs = [data[m, s, c] for s, c in itertools.product(group, ids)]
            full = [data["parent", s, c] for s, c in itertools.product(group, ids)]
            peer = [
                data["short_parent", s, c] for s, c in itertools.product(group, ids)
            ]
            depth = [
                int(x["hunt_evals"]) + int(x["second_evals"]) - int(x["cap"])
                for x in runs
            ]
            totals = {
                k: sum(int(r[k]) for r in runs)
                for k in (
                    "elapsed_ns",
                    "joint",
                    *COUNT,
                    "explicit_inputs",
                    "deferred_replaced",
                    "modes",
                )
            }
            costs.append(
                {
                    "receiver": rec,
                    "method": m,
                    "timing_block": 2 if m.startswith("full_") else 1,
                    "cells": len(runs),
                    "wall_s": round(totals.pop("elapsed_ns") / 1e9, 9),
                    **totals,
                    "wall_no_more_than_parent": sum(
                        int(l["elapsed_ns"]) <= int(r["elapsed_ns"])
                        for l, r in zip(runs, full, strict=True)
                    ),
                    "joint_no_more_than_parent": sum(
                        int(l["joint"]) <= int(r["joint"])
                        for l, r in zip(runs, full, strict=True)
                    ),
                    "hunt_plus_second_at_or_over_cap": sum(x >= 0 for x in depth),
                    "min_cap_excess": min(depth),
                    "max_cap_excess": max(depth),
                    "first_xyz_changed_vs_parent": sum(
                        l["first_xyz_sha256"] != r["first_xyz_sha256"]
                        for l, r in zip(runs, full, strict=True)
                    ),
                    "first_xyz_changed_vs_short_parent": sum(
                        l["first_xyz_sha256"] != r["first_xyz_sha256"]
                        for l, r in zip(runs, peer, strict=True)
                    ),
                }
            )
        for fold, subset in itertools.product((3, 10), ("all", "equal_heavy_charge")):
            edges = [
                (i, j)
                for i, j in itertools.permutations(range(n), 2)
                if k[j] >= fold * k[i] and (subset == "all" or desc[i] == desc[j])
            ]
            strong, weak = np.asarray(edges).T
            credit = {m: 1 + np.sign(s[weak] - s[strong]) for m, s in scores.items()}
            for m, s in scores.items():
                metrics.append(
                    {
                        "receiver": rec,
                        "fold": fold,
                        "subset": subset,
                        "method": m,
                        "n_ligands": n,
                        "n_pairs": len(edges),
                        "credit": float(credit[m].sum()) / 2,
                        "score_ties": int((credit[m] == 1).sum()),
                        "pKd_rho_neg_score": round(
                            float(spearmanr(pk, -s).statistic), 8
                        ),
                        "pKd_tau_b_neg_score": round(
                            float(kendalltau(pk, -s).statistic), 8
                        ),
                    }
                )
            if fold == 3:
                assert len(edges) == (403 if subset == "all" else 35)
                for a, b in CONTRASTS:
                    delta = credit[b] - credit[a]
                    for i in range(n):
                        retain = (strong != i) & (weak != i)
                        deletions.append(
                            {
                                "receiver": rec,
                                "subset": subset,
                                "reference": a,
                                "test": b,
                                "omitted": ids[i],
                                "n_pairs": int(retain.sum()),
                                "reference_credit": float(credit[a][retain].sum() / 2),
                                "test_credit": float(credit[b][retain].sum() / 2),
                                "delta_credit": float(delta[retain].sum() / 2),
                            }
                        )
                    if subset != "all":
                        continue
                    for h in np.flatnonzero(delta):
                        i, j = edges[h]
                        flips.append(
                            {
                                "receiver": rec,
                                "reference": a,
                                "test": b,
                                "strong": ids[i],
                                "weak": ids[j],
                                "Kd_strong_uM": k[i],
                                "Kd_weak_uM": k[j],
                                "same_size_charge": desc[i] == desc[j],
                                "same_scaffold": manifest[i]["scaffold"]
                                == manifest[j]["scaffold"],
                                "from_credit": credit[a][h] / 2,
                                "to_credit": credit[b][h] / 2,
                            }
                        )
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in (
        ("cells", list(data.values())),
        ("predictions", estimates),
        ("cost", costs),
        ("metrics", metrics),
        ("flips", flips),
        ("leave_one_out", deletions),
    ):
        write(out / (name + ".csv"), rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    analyze(a.work, a.out)
