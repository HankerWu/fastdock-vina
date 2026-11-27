"""Conditional paired uncertainty for the separate five-compound FABP4 assay."""

import csv
import itertools
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
with (HERE.parent / "results/assay_predictions.csv").open() as stream:
    scores = list(csv.DictReader(stream))
rows = []
for receptor in sorted({r["receiver"] for r in scores}):
    panel = {(r["arm"], r["id"]): r for r in scores if r["receiver"] == receptor}
    ligands = sorted({name for arm, name in panel})
    n = len(ligands)
    # Only five compounds: enumerate every equiprobable ordered resample.
    sample = np.array(list(itertools.product(range(n), repeat=n)))

    def centered(observed):
        ranks = rankdata(np.array(observed)[sample], axis=1)
        return ranks - ranks.mean(axis=1, keepdims=True)

    potency = centered([float(panel["old", name]["pK"]) for name in ligands])
    correlations = {}
    for arm in ["old", "new"]:
        predictions = centered(
            [-float(panel[arm, name]["top1_mean"]) for name in ligands]
        )
        denom = np.linalg.norm(potency, axis=1) * np.linalg.norm(predictions, axis=1)
        correlations[arm] = np.divide(
            (potency * predictions).sum(axis=1),
            denom,
            out=np.full(len(sample), np.nan),
            where=denom > 0,
        )
    delta = correlations["new"] - correlations["old"]
    lo, hi = np.nanpercentile(delta, [2.5, 97.5])
    rows.append(
        dict(
            receiver=receptor,
            ligands=n,
            resamples=len(sample),
            valid=np.isfinite(delta).sum(),
            delta_rho_p025=lo,
            delta_rho_p975=hi,
        )
    )
with (HERE / "results/fabp4_ligand_bootstrap.csv").open("w", newline="") as stream:
    w = csv.DictWriter(stream, rows[0], lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
