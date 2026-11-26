"""Ligand-only negative controls. Use the same 37 unique, 3-fold-separated pairs."""

import csv
import math
import statistics
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors

from run_panel import CORE, SERIES

dl = {}
for i, (ring, *_) in SERIES.items():
    mol = Chem.MolFromSmiles(CORE + ring)
    dl[i] = (
        Descriptors.MolLogP(mol),
        Descriptors.HeavyAtomCount(mol),
        Descriptors.MolWt(mol),
    )
with open("pairs.csv") as f:
    pairs = [
        r for r in csv.DictReader(f) if r["subset"] == "all" and r["score"] == "native"
    ]
with open("cells.csv") as f:
    cells = list(csv.DictReader(f))
gaps = {}
for score in ["native", "exclusion"]:
    e = {
        (r["ligand"], r["rec"], r["seed"]): float(r["energy"])
        for r in cells
        if r["score"] == score
    }
    for r in pairs:
        a, b = r["strong"], r["weak"]
        gaps[score, a, b] = statistics.median(
            v - e[a, rec, s] for (lig, rec, s), v in e.items() if lig == b
        )
rows = []
for label in (
    "native",
    "exclusion",
    "cLogP",
    "heavy atoms",
    "mass",
    "native plus octanol-transfer surrogate",
):
    for r in pairs:
        a, b = r["strong"], r["weak"]
        if label in ("native", "exclusion"):
            g = gaps[label, a, b]
        elif label == "native plus octanol-transfer surrogate":
            g = gaps["native", a, b] + (
                dl[a][0] - dl[b][0]
            ) * 0.00198720425864 * 298.15 * math.log(10)
        else:
            k = {"cLogP": 0, "heavy atoms": 1, "mass": 2}[label]
            g = dl[a][k] - dl[b][k]
        rows.append(
            {
                "method": label,
                "strong": a,
                "weak": b,
                "direction": g,
                "credit": 1 if g > 1e-9 else 0.5 if abs(g) < 1e-9 else 0,
            }
        )
with Path("controls.csv").open("w") as h:
    w = csv.DictWriter(h, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)
for meth in sorted({r["method"] for r in rows}):
    print(
        meth,
        [
            (
                sum(
                    r["credit"]
                    for r in rows
                    if r["method"] == meth
                    and (not same_charge or "4" not in (r["weak"], r["strong"]))
                ),
                sum(
                    1
                    for r in rows
                    if r["method"] == meth
                    and (not same_charge or "4" not in (r["weak"], r["strong"]))
                ),
            )
            for same_charge in (False, True)
        ],
    )
