"""Compound-level ranking, paired uncertainty and direction changes from fixed coordinates."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

RT = 0.00198720425864083 * 298.15


def summarize(frame, docks, out, detail_nodes=1024):
    rows = []
    edges = []
    ligand_rows = []
    for (count, target), p in frame.groupby(["nodes", "target"], sort=True):
        p = p.copy()
        p["separate"] = p["rigid"] + p["hinge"]
        arms = ["joint", "rigid", "hinge", "separate", "no_free", "min_delta"]
        for arm in arms:
            p[arm] += p["native"]
        p["size"] = -p["n_heavy"]
        arms = ["native", "size"] + arms
        control = docks[docks.target.eq(target)] if docks is not None else None
        if control is not None and not control.empty:
            wide = control.pivot(
                index=["ligand", "state", "seed"], columns="arm", values="score"
            ).reset_index()
            p = p.merge(wide, on=["ligand", "state", "seed"], validate="one_to_one")
            np.testing.assert_array_equal(p.pro, p.native)
            arms += ["pro_double", "upstream"]
        scores = (
            p.groupby(["ligand", "state"])[arms + ["kd_nM"]]
            .median()
            .groupby("ligand")
            .mean()
            .sort_values("kd_nM")
        )
        strong, weak = np.triu_indices(len(scores), 1)
        kd = scores.kd_nM.to_numpy()
        eligible = kd[weak] >= 3 * kd[strong]
        strong, weak = strong[eligible], weak[eligible]
        if count == detail_nodes:
            ligand_rows.extend(
                scores.assign(target=target).reset_index().to_dict("records")
            )
        exp = RT * np.log(kd[weak] / kd[strong])
        look = scores.index.to_numpy()
        rng = np.random.default_rng(29501)
        weights = rng.multinomial(
            len(scores), np.full(len(scores), 1 / len(scores)), size=5000
        )
        ew = weights[:, strong] * weights[:, weak]
        ew = ew[ew.sum(axis=1) > 0]
        for arm in arms:
            pred = scores[arm].to_numpy()
            gap = pred[weak] - pred[strong]
            credit = (gap > 0).astype(float) + 0.5 * (gap == 0)
            if arm == "native":
                base = credit.copy()
            delta = credit - base
            change = (ew @ delta) / ew.sum(axis=1)
            row = dict(
                nodes=count,
                target=target,
                arm=arm,
                ligands=len(scores),
                pairs=len(credit),
                credit=credit.sum(),
                gained=np.maximum(delta, 0).sum(),
                lost=-np.minimum(delta, 0).sum(),
                rho=spearmanr(pred, kd).statistic,
                tau_b=kendalltau(pred, kd).statistic,
                gap_mae=np.nan if arm == "size" else np.abs(gap - exp).mean(),
                bootstrap_delta_lo=np.quantile(change, 0.025),
                bootstrap_delta_hi=np.quantile(change, 0.975),
            )
            rows.append(row)
            if count == detail_nodes:
                for a, b, ex, pr, cr in zip(look[strong], look[weak], exp, gap, credit):
                    edges.append(
                        dict(
                            target=target,
                            arm=arm,
                            strong=a,
                            weak=b,
                            experimental_gap=ex,
                            predicted_gap=pr,
                            credit=cr,
                        )
                    )
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (("summary", rows), ("pairs", edges), ("ligands", ligand_rows)):
        pd.DataFrame(data).to_csv(
            out / (name + ".csv"), index=False, float_format="%.9g"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("cells", nargs="+")
    parser.add_argument("--docks", type=Path)
    parser.add_argument("--detail-nodes", type=int, default=1024)
    a = parser.parse_args()
    summarize(
        pd.concat([pd.read_csv(path) for path in a.cells]),
        pd.read_csv(a.docks) if a.docks else None,
        a.out,
        a.detail_nodes,
    )
