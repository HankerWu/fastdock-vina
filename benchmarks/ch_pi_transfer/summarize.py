"""Recompute paired geometric and within-assay results from the mode-level output."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
ARMS = ("native", "ligand_CH", "protein_CH", "bidirectional")
AFFINITY_ARMS = (*ARMS, "rescore_native_winner", "native_at_changed_winner")
KEY = ["target", "pdb", "compound", "cap", "seed"]


def describe(x):
    q = np.quantile(x, [0.25, 0.5, 0.75, 0.90, 0.95])
    return {
        "mean": np.mean(x),
        "q25": q[0],
        "median": q[1],
        "q75": q[2],
        "p90": q[3],
        "p95": q[4],
        "maximum": max(x),
    }


def geometry(poses):
    cells = []
    for key, pool in poses.groupby(KEY, sort=True):
        for arm in ARMS:
            selected = pool.sort_values([arm, "mode"])
            top = selected.iloc[0]
            oracle = pool.rmsd.min()
            cells.append(
                dict(
                    zip(KEY, key),
                    arm=arm,
                    modes=len(pool),
                    top_mode=top["mode"],
                    score=top[arm],
                    top1=top.rmsd,
                    best5=selected.iloc[:5].rmsd.min(),
                    pool_min=oracle,
                    regret=top.rmsd - oracle,
                )
            )
    c = pd.DataFrame(cells)
    base = c[c.arm == "native"][KEY + ["top1", "top_mode"]].rename(
        columns={"top1": "native_top1", "top_mode": "native_mode"}
    )
    c = c.merge(base, on=KEY, validate="many_to_one")
    c["delta"] = c.top1 - c.native_top1
    rows = []
    for scope, by in [("target", ["target"]), ("structure", ["target", "pdb"])]:
        for keys, g in c.groupby(by + ["arm"], sort=True):
            ids = dict(zip(by + ["arm"], keys))
            for metric in ["top1", "best5", "pool_min", "regret", "delta"]:
                rows.append(
                    dict(
                        **ids,
                        scope=scope,
                        measure=metric,
                        n=len(g),
                        changed=(g.native_mode != g.top_mode).sum(),
                        improved=(g.delta < 0).sum(),
                        worse=(g.delta > 0).sum(),
                        **describe(g[metric]),
                    )
                )
    return c, pd.DataFrame(rows)


def ranks(kd, err, scores):
    """Vectorized per-sample Spearman and tau-b with proper ties; pair cutoff in Kd space."""
    kd = np.atleast_2d(kd)
    scores = np.atleast_2d(scores)
    err = np.atleast_2d(err)
    i, j = np.triu_indices(kd.shape[1], k=1)
    # Score and Kd are both smaller for a stronger ligand. pKd and -score
    # give the same rank statistic. Raw Kd is needed for the gap/uncertainty.
    k, s = np.sign(kd[:, i] - kd[:, j]), np.sign(scores[:, i] - scores[:, j])
    direction = k * s
    different = k != 0
    triple = (np.maximum(kd[:, i], kd[:, j]) >= 3 * np.minimum(kd[:, i], kd[:, j])) & (
        (kd[:, i] - err[:, i] > kd[:, j] + err[:, j])
        | (kd[:, j] - err[:, j] > kd[:, i] + err[:, i])
    )
    xr, yr = rankdata(kd, axis=1), rankdata(scores, axis=1)
    x, y = xr - xr.mean(axis=1, keepdims=True), yr - yr.mean(axis=1, keepdims=True)

    def credit(mask):
        return (((direction > 0) + 0.5 * (s == 0)) * mask).sum(axis=1)

    with np.errstate(invalid="ignore", divide="ignore"):
        return {
            "rho": (x * y).sum(1) / np.sqrt((x * x).sum(1) * (y * y).sum(1)),
            "tau": direction.sum(1) / np.sqrt(different.sum(1) * (s != 0).sum(1)),
            "pairs": different.sum(1),
            "concordance": credit(different) / different.sum(1),
            "resolved": triple.sum(1),
            "resolved_credit": credit(triple),
            "resolved_concordance": credit(triple) / triple.sum(1),
        }


def affinity(c, poses):
    # Fixed-winner controls separate reordering from adding a ligand-specific offset.
    control = []
    for key, pool in poses.groupby(KEY, sort=True):
        native = pool.loc[pool.native.idxmin()]
        changed = pool.loc[pool.bidirectional.idxmin()]
        for arm, score in [
            ("rescore_native_winner", native.bidirectional),
            ("native_at_changed_winner", changed.native),
        ]:
            control.append(dict(zip(KEY, key), arm=arm, score=score))
    c = pd.concat([c, pd.DataFrame(control)], ignore_index=True)
    label = pd.read_csv(HERE.parent / "galectin_affinity/table1.csv").set_index(
        "compound"
    )
    target_rows, ligand_rows, intervals, pair_rows = [], [], [], []
    for target in ["Gal1", "Gal3"]:
        series = c[c.target == target].copy()
        series["compound"] = series.compound.astype(int)
        med = series.groupby(["pdb", "compound", "arm"]).score.median().unstack()
        by_state = {p: med.xs(p, level="pdb") for p in series.pdb.unique()}
        by_state["mean_of_states"] = med.groupby(level="compound").mean()
        for pdb, table in by_state.items():
            # Means of medians of millesimal Vina reports can be exactly tied
            # in decimal arithmetic. Discard binary round-off below 1e-9.
            table = table.round(9)
            kd = label.loc[table.index, target + "_Kd_uM"].values
            err = label.loc[table.index, target + "_error_uM"].values
            records = table.assign(
                kd_uM=kd, error_uM=err, pKd=6 - np.log10(kd), target=target, pdb=pdb
            ).reset_index()
            ligand_rows.extend(records.to_dict("records"))
            rng = np.random.default_rng(271828)
            draws = rng.integers(0, len(kd), (4000, len(kd)))
            refstats = ranks(kd[draws], err[draws], table.native.values[draws])
            for arm in AFFINITY_ARMS:
                score = table[arm].values
                measure = ranks(kd, err, score)
                target_rows.append(
                    dict(
                        target=target,
                        pdb=pdb,
                        arm=arm,
                        n=len(kd),
                        **{k: v[0] for k, v in measure.items()},
                    )
                )
                if arm != "bidirectional":
                    continue
                boot = ranks(kd[draws], err[draws], score[draws])
                for key in ["rho", "tau", "concordance", "resolved_concordance"]:
                    d = boot[key] - refstats[key]
                    valid = d[np.isfinite(d)]
                    intervals.append(
                        {
                            "target": target,
                            "pdb": pdb,
                            "measure": key,
                            "draws": len(d),
                            "valid": len(valid),
                            "lo": np.quantile(valid, 0.025),
                            "hi": np.quantile(valid, 0.975),
                        }
                    )
            for i, j in zip(*np.triu_indices(len(kd), k=1)):
                lo, hi = sorted([i, j], key=kd.__getitem__)
                pair_rows.append(
                    dict(
                        target=target,
                        pdb=pdb,
                        lig1=table.index[lo],
                        lig2=table.index[hi],
                        kd1=kd[lo],
                        kd2=kd[hi],
                        experimental_tie=kd[hi] == kd[lo],
                        resolved=kd[hi] >= 3 * kd[lo]
                        and kd[hi] - err[hi] > kd[lo] + err[lo],
                        **{
                            arm + "_weak_minus_strong": table[arm].iloc[hi]
                            - table[arm].iloc[lo]
                            for arm in ARMS
                        },
                    )
                )
    return {
        "affinity": target_rows,
        "ligand_scores": ligand_rows,
        "affinity_intervals": intervals,
        "pairs": pair_rows,
    }


def run(output):
    poses = pd.read_csv(output / "poses.csv")
    cells, geom = geometry(poses)
    tables = dict(cells=cells, geometry=geom, **affinity(cells, poses))
    # Same compound at all seeds and both receptor conformations travels together.
    g = cells[(cells.arm == "bidirectional") & cells.target.isin(["Gal1", "Gal3"])]
    cis = []
    for target, panel in g.groupby("target"):
        per_ligand = panel.groupby("compound").delta.mean().values
        rng = np.random.default_rng(271828)
        bootstrap = rng.choice(per_ligand, (4000, len(per_ligand))).mean(axis=1)
        cis.append(
            {
                "target": target,
                "n_compounds": len(per_ligand),
                "mean": per_ligand.mean(),
                "low": np.quantile(bootstrap, 0.025),
                "high": np.quantile(bootstrap, 0.975),
            }
        )
    tables["geometry_intervals"] = cis
    # Each distinct error is an ECDF jump. Repeated values retain their full mass.
    ecdf = []
    for (target, arm), panel in cells.groupby(["target", "arm"]):
        for measure in ["top1", "best5", "pool_min", "regret", "delta"]:
            x, counts = np.unique(panel[measure], return_counts=True)
            for a, p in zip(x, counts.cumsum() / counts.sum()):
                ecdf.append(
                    {
                        "target": target,
                        "arm": arm,
                        "measure": measure,
                        "error_A": a,
                        "fraction_at_or_below": p,
                    }
                )
    tables["ecdf"] = ecdf
    for name, table in tables.items():
        pd.DataFrame(table).to_csv(
            output / f"{name}.csv", index=False, float_format="%.9f"
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    run(p.parse_args().out)
