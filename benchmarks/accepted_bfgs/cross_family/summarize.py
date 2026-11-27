"""Aggregate all prespecified matched starts, retaining per-pair and per-target comparisons."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from replay import STARTS, TARGETS, coords, write_rows


def describe(values):
    x = np.array(values)
    quantiles = np.quantile(x, [0.25, 0.5, 0.75, 0.9, 0.95])
    return {
        "n": len(x),
        "mean": float(x.mean()),
        **dict(
            zip(
                ("q25", "median", "q75", "p90", "p95"),
                map(float, quantiles),
                strict=True,
            )
        ),
        "max": float(x.max()),
    }


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def report(work):
    metrics = (
        "end_rmsd_A",
        "end_zero_score",
        "end_objective",
        "returned_native_objective",
        "native_fg",
    )
    results = read_rows(work / "rows.csv")
    observed = {
        (r["pdb"], r["start"], r["arm"]): r for r in read_rows(work / "native.csv")
    }
    by_key = {}
    for row in results:
        key = (row["pdb"], row["start"], row["arm"])
        trace = observed[key]
        if key in by_key or row["end_sha256"] != trace["verified_end_sha256"]:
            raise ValueError(f"Missing production parity: {key}")
        by_key[key] = {**row, **trace}
    if len(by_key) != len(observed) or len(by_key) != 48:
        raise ValueError("Incomplete paired panel")
    pairs = []
    for pdb in TARGETS:
        for start in STARTS:
            a, b = (by_key[pdb, start, arm] for arm in ("parent", "fixed"))
            for k in (
                "start_sha256",
                "start_rmsd_A",
                "start_zero_score",
                "start_objective",
            ):
                if a[k] != b[k]:
                    raise ValueError(f"Unmatched {k} for {pdb} {start}")
            pair = {
                k: a[k]
                for k in (
                    "pdb",
                    "target",
                    "start",
                    "start_rmsd_A",
                    "start_zero_score",
                    "start_objective",
                )
            }
            pair["identical_output"] = int(a["end_sha256"] == b["end_sha256"])
            for k in metrics:
                left, right = float(a[k]), float(b[k])
                pair.update(
                    {
                        f"parent_{k}": left,
                        f"fixed_{k}": right,
                        f"delta_{k}": right - left,
                    }
                )
            pairs.append(pair)
    write_rows(work / "paired.csv", pairs)
    groups = {}
    for group, subset in (
        ("all", pairs),
        ("projected", [p for p in pairs if p["start"] == "projected"]),
        ("wrong", [p for p in pairs if p["start"] != "projected"]),
    ):
        groups[group] = {}
        for arm in ("parent", "fixed"):
            values = [p[f"{arm}_end_rmsd_A"] for p in subset]
            groups[group][arm] = {
                **describe(values),
                "within_2A": sum(v <= 2 for v in values),
                "above_5A": sum(v > 5 for v in values),
                "native_fg": sum(int(p[f"{arm}_native_fg"]) for p in subset),
            }
        groups[group]["paired_delta_rmsd_A"] = describe(
            [p["delta_end_rmsd_A"] for p in subset]
        )
    rank = []
    for pdb in TARGETS:
        for arm in ("parent", "fixed"):
            projected = by_key[pdb, "projected", arm]
            for start in STARTS[1:]:
                wrong = by_key[pdb, start, arm]
                diff = float(wrong["end_zero_score"]) - float(
                    projected["end_zero_score"]
                )
                rank.append(
                    {
                        "pdb": pdb,
                        "target": projected["target"],
                        "wrong_start": start,
                        "arm": arm,
                        "wrong_minus_projected_zero_score": diff,
                    }
                )
    write_rows(work / "score_pairs.csv", rank)
    summary = {
        "rmsd_A": groups,
        "identical_output_files": sum(p["identical_output"] for p in pairs),
        "better_projected_score_out_of_16": {
            arm: sum(
                r["wrong_minus_projected_zero_score"] > 0
                for r in rank
                if r["arm"] == arm
            )
            for arm in ("parent", "fixed")
        },
        "timing": {
            arm: {
                k: sum(float(r[k]) for r in results if r["arm"] == arm)
                for k in ("cpu_s", "wall_s")
            }
            for arm in ("parent", "fixed")
        },
    }
    summary["unchanged_from_input_at_pdbqt_precision"] = {
        arm: sum(
            coords(work / "starts" / p["pdb"] / (p["start"] + ".pdbqt"))
            == coords(work / "raw" / f"{p['pdb']}-{p['start']}-{arm}.pdbqt")
            for p in pairs
        )
        for arm in ("parent", "fixed")
    }
    (work / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    report(ap.parse_args().work)
