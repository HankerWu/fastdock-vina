"""Reconcile two independent simulation tables with Morton 29 C ITC and paired native scores."""

import csv
import itertools
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

from rdkit import Chem

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
RT = 0.00198720425864 * 302.15
METHODS = (
    "cm2_single",
    "cm2_multi",
    "cm2_val",
    "bcc_val",
    "srti_bcc",
    "hrex_bcc",
    "srti_mixed",
    "hrex_mixed",
    "native",
)


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def canonical(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None or Chem.GetFormalCharge(m) != 0:
        raise ValueError(smiles)
    return Chem.MolToSmiles(m)


def value(row, method):
    if method.endswith("_mixed"):
        resp = method.replace("_mixed", "_resp")
        method = resp if row[resp] else method.replace("_mixed", "_bcc")
    return float(row[method]), float(row[method + "_sd"])


def reconcile():
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    assay = {r["pdb"]: r for r in read(BASE / "cutoff_switch/t4l/affinity_check.csv")}
    graphs = json.loads((BASE / "cutoff_switch/t4l/alignment.json").read_text())
    source = {r["pdb"]: r for r in read(HERE / "published.csv")}
    if source.keys() != assay.keys():
        raise ValueError("Incomplete panel")
    natives = defaultdict(dict)
    for r in read(BASE / "anchor_search/transfer_holdout_t4l/cells.csv"):
        if r["arm"] == "native":
            natives[r["ligand"]][r["receptor"], r["seed"]] = float(r["score"])
    for p, r in source.items():
        if r["ligand"] != assay[p]["name"] or canonical(r["smiles"]) != canonical(
            graphs[p]["smiles"]
        ):
            raise ValueError(f"Chemical identity mismatch: {p}")
        if abs(float(r["exp"]) - float(assay[p]["ITC_dG_kcal_mol"])) > 1e-8:
            raise ValueError(f"Reported experimental free energies disagree: {p}")
    pairs = []
    for a, b in itertools.combinations(
        sorted(assay, key=lambda x: float(assay[x]["Kd_uM"])), 2
    ):
        ratio = float(assay[b]["Kd_uM"]) / float(assay[a]["Kd_uM"])
        if len(natives[a]) != 4 or natives[a].keys() != natives[b].keys():
            raise ValueError("Unpaired native blocks")
        for arm in METHODS:
            if arm == "native":
                diff = statistics.median(
                    natives[b][k] - natives[a][k] for k in natives[b]
                )
                sd0 = sd1 = None
            else:
                x, s0 = value(source[a], arm)
                y, s1 = value(source[b], arm)
                diff = y - x
                # Published entries can share simulation legs. True pair SD is between these bounds.
                sd0 = abs(s0 - s1)
                sd1 = s0 + s1
            observed = RT * math.log(ratio)
            pairs.append(
                {
                    "strong": a,
                    "weak": b,
                    "kd_ratio": ratio,
                    "threefold": int(ratio >= 3),
                    "arm": arm,
                    "experimental_gap": observed,
                    "experimental_error_sum": float(assay[a]["reported_dG_error"])
                    + float(assay[b]["reported_dG_error"]),
                    "calculated_gap": diff,
                    "pair_error": diff - observed,
                    "credit": int(diff > 0) + 0.5 * (diff == 0),
                    "sd_min": sd0,
                    "sd_max": sd1,
                }
            )
    write(out / "published_pairs.csv", pairs)
    summaries = []
    for arm, restricted in itertools.product(METHODS, (False, True)):
        rows = [
            r for r in pairs if r["arm"] == arm and (not restricted or r["threefold"])
        ]
        robust = [
            r
            for r in rows
            if r["sd_max"] is not None and abs(r["calculated_gap"]) > r["sd_max"]
        ]
        loo = []
        for pdb in source:
            remain = [r["credit"] for r in rows if pdb not in (r["strong"], r["weak"])]
            loo.append(statistics.mean(remain))
        summaries.append(
            {
                "arm": arm,
                "threshold": "3fold" if restricted else "all",
                "n": len(rows),
                "credit": sum(r["credit"] for r in rows),
                "gap_mae": statistics.mean(abs(r["pair_error"]) for r in rows),
                "gap_rmse": math.sqrt(
                    statistics.mean(r["pair_error"] ** 2 for r in rows)
                ),
                "robust_n": len(robust) if arm != "native" else None,
                "robust_credit": sum(r["credit"] for r in robust)
                if arm != "native"
                else None,
                "omit_one_ligand_min": min(loo),
                "omit_one_ligand_max": max(loo),
            }
        )
    write(out / "published_summary.csv", summaries)
    apo = []
    for p, r in source.items():
        if r["holo_cm2_multi"]:
            apo.append(
                {
                    "pdb": p,
                    "ligand": r["ligand"],
                    "apo_multi": r["cm2_multi"],
                    "holo_multi": r["holo_cm2_multi"],
                    "apo_released": r["cm2_val"],
                    "holo_released": r["holo_cm2_val"],
                    "initial_holo_minus_apo": float(r["holo_cm2_multi"])
                    - float(r["cm2_multi"]),
                    "released_holo_minus_apo": float(r["holo_cm2_val"])
                    - float(r["cm2_val"]),
                }
            )
    write(out / "apo_holo.csv", apo)


if __name__ == "__main__":
    reconcile()
