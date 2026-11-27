"""Propagate the finite free-library minimum through unchanged pose and assay scores."""

import argparse
import csv
import itertools
import lzma
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

SOURCE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE.parent / "strain_components"))
from summarize import pair_rows, write

ARMS = {"full": "vacuum", "gb_reselected": "solution"}


def read(path):
    with lzma.open(path, "rt") if path.suffix == ".xz" else path.open() as stream:
        return list(csv.DictReader(stream))


def finite_references(path):
    groups = defaultdict(list)
    for row in read(path):
        groups[row["compound"]].append(row)
    return groups


def offsets(refs, new_refs, scales, extra):
    rows, shifts = [], {}
    if refs.keys() != new_refs.keys() or refs.keys() != scales.keys():
        raise ValueError("Ligand panels disagree")
    for c, library in sorted(refs.items()):
        if len(library) != 8:
            raise ValueError("The source ensemble must have eight converged samples")
        for arm, field in ARMS.items():
            energies = sorted(float(r[field]) for r in library)
            f8 = min(energies)
            g = scales[c]
            shifts4 = np.array(
                [
                    g * (f8 - min(subset))
                    for subset in itertools.combinations(energies, 4)
                ]
            )
            fresh = min(float(r[field]) for r in new_refs[c])
            shift = g * (f8 - fresh)
            shift32 = None
            if c in extra and field == "vacuum":
                if abs(float(extra[c]["free_8"]) - f8) > 1e-6:
                    raise ValueError("Incompatible historical vacuum references")
                shift32 = g * (f8 - float(extra[c]["free_32"]))
            rows.append(
                {
                    "compound": c,
                    "arm": arm,
                    "scale": g,
                    "historical_energy": f8,
                    "min4_energy_max": energies[-4],
                    "fresh_energy": fresh,
                    "fresh_converged": len(new_refs[c]),
                    "fresh_score_shift": shift,
                    "score_shift_32": shift32,
                }
            )
            shifts[c, arm] = (shifts4, shift, shift32)
    return rows, shifts


def fixed_pairs(modes, meta):
    # Keep the original panels' distinct receptor/seed aggregation definitions.
    groups = defaultdict(lambda: defaultdict(dict))
    for row in modes:
        if row["assay"] == "lectin":
            continue
        assay, rec, c = (row[k] for k in ("assay", "pdb", "compound"))
        for target in (rec, "both") if assay == "CHEMBL661804" else (rec,):
            groups[assay, target][c][rec, row["seed"]] = row
    pairs = []
    for (assay, rec), states in sorted(groups.items()):
        for r in pair_rows(assay, rec, states, meta, ("native", *ARMS)):
            a, b = r["stronger"], r["weaker"]
            matched = (
                meta[a]["charge"] == meta[b]["charge"]
                and meta[a]["heavy"] == meta[b]["heavy"]
            )
            pairs.append(dict(**r, size_charge_matched=matched))
    lectin = defaultdict(dict)
    for r in read(SOURCE / "pairs.csv.xz"):
        if r["arm"] in ("native", *ARMS):
            key = (r["target"], r["pdb"], r["stronger"], r["weaker"], r["fold"])
            lectin[key][r["arm"] + "_gap"] = float(r["gap"])
            lectin[key]["size_charge_matched"] = r["same_formula"] == "True"
    for (assay, rec, a, b, fold), row in lectin.items():
        pairs.append(
            dict(assay=assay, receptor=rec, stronger=a, weaker=b, fold=fold, **row)
        )
    return pairs


def credit(gap, tol):
    return 0.5 if abs(gap) < tol else float(gap > 0)


def evaluate(pairs, shifts, independent32):
    rows = []
    for p in pairs:
        a, b = p["stronger"], p["weaker"]
        tol = 1e-9 if p["assay"] in ("Gal1", "Gal3") else 1e-10
        native = credit(p["native_gap"], tol)
        for arm in ARMS:
            sa, sb = shifts[a, arm], shifts[b, arm]
            gap = p[arm + "_gap"]
            historical = credit(gap, tol)
            # Positive gaps favor the experimentally stronger compound.
            gaps4 = gap + sb[0][:, None] - sa[0][None, :]
            credits4 = np.where(abs(gaps4) < tol, 0.5, gaps4 > 0)
            fresh_gap = gap + sb[1] - sa[1]
            fresh = credit(fresh_gap, tol)
            key32 = (p["assay"], p["receptor"], a, b)
            new32 = None
            if arm == "full" and key32 in independent32:
                gap32 = gap + sb[2] - sa[2]
                new32 = credit(gap32, tol)
                if abs(gap32 - independent32[key32]) > 1e-7:
                    raise ValueError("Independent 32-conformer score does not match")
            rows.append(
                {
                    "assay": p["assay"],
                    "receptor": p["receptor"],
                    "stronger": a,
                    "weaker": b,
                    "fold": p["fold"],
                    "matched": p["size_charge_matched"],
                    "arm": arm,
                    "native": native,
                    "historical_gap": gap,
                    "historical": historical,
                    "min_gap4": gaps4.min(),
                    "max_gap4": gaps4.max(),
                    "expected_credit4": credits4.mean(),
                    "invariant": bool(np.all(credits4 == historical)),
                    "fresh_gap": fresh_gap,
                    "fresh": fresh,
                    "credit32": new32,
                }
            )
    return rows


def summarize(rows):
    groups = defaultdict(list)
    for r in rows:
        for subset in ("all", "matched") if r["matched"] else ("all",):
            groups[r["assay"], r["receptor"], r["arm"], subset].append(r)
    summary = []
    for (assay, rec, arm, subset), rr in sorted(groups.items()):
        stable = [r for r in rr if r["invariant"]]
        fragile = [r for r in rr if not r["invariant"]]

        def changed(group, method):
            return sum(r[method] != r["historical"] for r in group)

        ligands = {r[k] for r in rr for k in ("stronger", "weaker")}
        deletions = [
            sum(
                r["fresh"] - r["native"]
                for r in rr
                if c not in (r["stronger"], r["weaker"])
            )
            for c in ligands
        ]
        summary.append(
            {
                "assay": assay,
                "receptor": rec,
                "arm": arm,
                "subset": subset,
                "pairs": len(rr),
                "native": sum(r["native"] for r in rr),
                "historical": sum(r["historical"] for r in rr),
                "expect4": round(sum(r["expected_credit4"] for r in rr), 4),
                "fresh": sum(r["fresh"] for r in rr),
                "gained": sum(r["fresh"] > r["historical"] for r in rr),
                "lost": sum(r["fresh"] < r["historical"] for r in rr),
                "invariant4": len(stable),
                "changed_fresh_invariant": changed(stable, "fresh"),
                "variable4": len(fragile),
                "changed_fresh_variable": changed(fragile, "fresh"),
                "credit32": sum(r["credit32"] for r in rr)
                if rr[0]["credit32"] is not None
                else None,
                "changed32_invariant": changed(stable, "credit32")
                if rr[0]["credit32"] is not None
                else None,
                "changed32_variable": changed(fragile, "credit32")
                if rr[0]["credit32"] is not None
                else None,
                "min_loo_net_credit_vs_native": min(deletions),
                "max_loo_net_credit_vs_native": max(deletions),
            }
        )
    return summary


def check_totals(summary):
    # Verify the full old database against separate archived publications.
    expected = {}
    for r in read(SOURCE / "corrected_summary.csv"):
        if r["subset"] == "all":
            for arm in ARMS:
                expected["CHEMBL661804", r["receptor"], arm] = (
                    float(r["native"]),
                    float(r[arm]),
                )
    for r in read(SOURCE / "protease_summary.csv"):
        if r["subset"] == "all" and r["assay"] != "CHEMBL661804":
            for arm in ARMS:
                expected[r["assay"], r["receptor"], arm] = (
                    float(r["native"]),
                    float(r[arm]),
                )
    lectin = {
        (r["target"], r["pdb"], r["arm"]): float(r["credit"])
        for r in read(SOURCE / "affinity.csv")
        if r["subset"] == "all"
    }
    for (target, rec, arm), score in lectin.items():
        if arm in ARMS:
            expected[target, rec, arm] = lectin[target, rec, "native"], score
    reconstructed = {
        (r["assay"], r["receptor"], r["arm"]): (r["native"], r["historical"])
        for r in summary
        if r["subset"] == "all"
    }
    if reconstructed != expected:
        raise ValueError("Historical scores or pair aggregation do not reproduce")


def main(out):
    modes = read(SOURCE / "modes.csv.xz")
    modes = [r for r in modes if r["compound"] != "CHEMBL369480"]
    modes.extend(read(SOURCE / "corrected_modes.csv"))
    scales = {}
    for r in modes:
        c, scale = r["compound"], float(r["scale"])
        if scales.setdefault(c, scale) != scale:
            raise ValueError("Reference offset does not commute with pose aggregation")
    extra = {
        r["compound"]: r for r in read(SOURCE.parent / "strain/free_references.csv")
    }
    meta = {
        r["compound"]: r for r in read(SOURCE.parent / "strain_components/scores.csv")
    }
    independent32 = {
        (r["assay"], r["receptor"], r["stronger"], r["weaker"]): float(r["gap"])
        for r in read(SOURCE.parent / "strain/fresh_mmff_pairs.csv.xz")
        if r["target"] == "BACE1" and r["arm"] == "unit_32" and r["policy"] == "first"
    }
    ligand_rows, shifts = offsets(
        finite_references(SOURCE / "references.csv.xz"),
        finite_references(out / "independent_references.csv.xz"),
        scales,
        extra,
    )
    rows = evaluate(fixed_pairs(modes, meta), shifts, independent32)
    summaries = summarize(rows)
    check_totals(summaries)
    write(out / "ligands.csv", ligand_rows)
    write(out / "pairs.csv.xz", rows)
    write(out / "summary.csv", summaries)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    main(parser.parse_args().out)
