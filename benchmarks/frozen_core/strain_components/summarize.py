"""Reconstruct assay-paired accuracy directly from compound scores and Ki measurements."""

import argparse
import csv
import itertools
import lzma
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from statistics import median

from scipy.stats import kendalltau, spearmanr

ARMS = ("native", "full", "without_elec", "elec_only")
HERE = Path(__file__).resolve().parent


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "wt", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def assignments(a, b):
    yield "all"
    if a["charge"] != b["charge"]:
        yield "cross_charge"
    else:
        yield "same_charge"
        yield "both_charge_" + a["charge"]
        if a["heavy"] == b["heavy"]:
            yield "same_heavy_charge"
    if a["scaffold"] == b["scaffold"]:
        yield "same_scaffold"
    else:
        yield "different_scaffold"


def pair_rows(assay, receptor, states, meta):
    for a, b in itertools.combinations(
        sorted(states, key=lambda c: Decimal(meta[c]["ki_nm"])), 2
    ):
        fold = Decimal(meta[b]["ki_nm"]) / Decimal(meta[a]["ki_nm"])
        if fold < 3:
            continue
        cells = sorted(states[a])
        if cells != sorted(states[b]):
            raise ValueError("Inconsistent receptor/seed blocks")
        rows = {
            "assay": assay,
            "receptor": receptor,
            "stronger": a,
            "weaker": b,
            "fold": float(fold),
        }
        for arm in ARMS:
            # BACE uses a difference of seed-medians, Factor Xa a median of
            # matched differences. Preserve both published definitions exactly.
            if assay == "CHEMBL661804":
                gap = median(
                    float(states[b][k][arm]) - float(states[a][k][arm]) for k in cells
                )
            else:
                gap = median(float(r[arm]) for r in states[b].values()) - median(
                    float(r[arm]) for r in states[a].values()
                )
            rows[arm + "_gap"] = gap
            rows[arm] = 0.5 if abs(gap) < 1e-10 else int(gap > 0)
        yield rows


def sensitivity(edges, ligands, arm, base):
    values = []
    for c in ligands:
        subset = [r for r in edges if c not in (r["stronger"], r["weaker"])]
        if subset:
            values.append(100 * sum(r[arm] - r[base] for r in subset) / len(subset))
    return round(min(values), 3), round(max(values), 3)


def summarize(out):
    by = defaultdict(lambda: defaultdict(dict))
    meta = {}
    for row in read(out / "scores.csv"):
        assay, c, rec, seed = (
            row[k] for k in ("assay", "compound", "receptor", "seed")
        )
        meta[c] = row
        for receiver in [rec, "both"] if assay == "CHEMBL661804" else [rec]:
            if (rec, seed) in by[assay, receiver][c]:
                raise ValueError("Duplicate cell")
            by[assay, receiver][c][rec, seed] = row

    pairs, summaries, correlations, influence = [], [], [], []
    for (assay, receiver), states in sorted(by.items()):
        rows = list(pair_rows(assay, receiver, states, meta))
        pairs.extend(rows)
        subsets = defaultdict(list)
        for row in rows:
            a, b = row["stronger"], row["weaker"]
            for subset in assignments(meta[a], meta[b]):
                subsets[subset].append(row)
            if row["fold"] >= 10:
                subsets["tenfold"].append(row)
        for subset, edges in sorted(subsets.items()):
            summaries.append(
                dict(
                    assay=assay,
                    receptor=receiver,
                    subset=subset,
                    pairs=len(edges),
                    **{arm: sum(r[arm] for r in edges) for arm in ARMS},
                )
            )
        for arm, base in (
            ("full", "native"),
            ("without_elec", "full"),
            ("elec_only", "native"),
        ):
            lo, hi = sensitivity(rows, states, arm, base)
            influence.append(
                {
                    "assay": assay,
                    "receptor": receiver,
                    "comparison": f"{arm}-{base}",
                    "n_ligands": len(states),
                    "min_loo_delta_pp": lo,
                    "max_loo_delta_pp": hi,
                }
            )
        if receiver == "both":
            continue  # The median of four pair gaps has no defined scalar ligand ranking.
        ligands = sorted(states)
        measured = [float(meta[c]["ki_nm"]) for c in ligands]
        for arm in ARMS:
            pred = [median(float(r[arm]) for r in states[c].values()) for c in ligands]
            correlations.append(
                {
                    "assay": assay,
                    "receptor": receiver,
                    "arm": arm,
                    "ligands": len(ligands),
                    "spearman": round(spearmanr(measured, pred).statistic, 6),
                    "kendall_tau_b": round(kendalltau(measured, pred).statistic, 6),
                }
            )

    originals = {}
    strain = HERE.parent / "strain"
    for r in read(strain / "fresh_mmff_summary.csv"):
        if (
            r["target"] == "BACE1"
            and r["arm"] in ("native", "unit_8")
            and r["policy"] == "first"
        ):
            key = (
                r["assay"],
                r["receptor"],
                "native" if r["arm"] == "native" else "full",
            )
            originals[key] = (int(r["pairs"]), float(r["credit"]))
    for r in read(strain / "fxa_summary.csv"):
        if (
            r["method"] in ("native", "full_unit")
            and r["stratum"] == "all"
            and r["policy"] == "first"
        ):
            key = (
                "CHEMBL661804",
                r["receptor"],
                "native" if r["method"] == "native" else "full",
            )
            originals[key] = (int(r["pairs"]), float(r["credit"]))
    replay = {
        (r["assay"], r["receptor"], arm): (r["pairs"], r[arm])
        for r in summaries
        if r["subset"] == "all"
        for arm in ("native", "full")
    }
    if replay != originals:
        raise ValueError(f"Archived totals diverge: {replay} versus {originals}")
    write(out / "pairs.csv.xz", pairs)
    write(out / "summary.csv", summaries)
    write(out / "correlations.csv", correlations)
    write(out / "influence.csv", influence)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    summarize(parser.parse_args().out)
