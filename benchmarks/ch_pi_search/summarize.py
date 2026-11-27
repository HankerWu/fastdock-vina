"""Report each protein and each nominal budget separately, with paired contrasts."""

from argparse import ArgumentParser
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

IDS = ["panel", "query", "receptor", "seed", "cap"]
ARMS = {
    ("native", "native"): "baseline",
    ("field", "native"): "search_only",
    ("native", "field"): "rerank_only",
    ("field", "field"): "both",
}


def summarize(directory):
    cells = pd.read_csv(directory / "cells.csv")
    cells["arm"] = [ARMS[k] for k in zip(cells.search, cells.rerank)]
    dist = []
    for key, rows in cells.groupby(["panel", "cap", "arm"]):
        for metric in ("top1", "pool_min", "regret", "best5"):
            x = rows[metric]
            dist.append(
                dict(
                    panel=key[0],
                    cap=key[1],
                    arm=key[2],
                    metric=metric,
                    n=len(x),
                    mean=x.mean(),
                    median=x.median(),
                    q25=x.quantile(0.25),
                    q75=x.quantile(0.75),
                    p90=x.quantile(0.90),
                    p95=x.quantile(0.95),
                    worst=x.max(),
                )
            )
    pairs = cells.pivot(index=IDS, columns="arm", values="top1").reset_index()
    for arm in ("baseline", "search_only", "rerank_only"):
        pairs["both_minus_" + arm] = pairs.both - pairs[arm]
    pairs["interaction"] = (
        pairs.both - pairs.search_only - pairs.rerank_only + pairs.baseline
    )
    contrasts = []
    for (panel, cap), rows in pairs.groupby(["panel", "cap"]):
        for col in [
            "both_minus_baseline",
            "both_minus_search_only",
            "both_minus_rerank_only",
            "interaction",
        ]:
            x = rows[col]
            # These are descriptive ranges for a small, fixed chemical series.
            blocks = rows.groupby("query")[col].mean()
            limits = (None, None)
            if panel in ("GBA1", "GAA"):
                draws = [
                    np.mean(sample) for sample in product(blocks, repeat=len(blocks))
                ]
                limits = np.quantile(draws, [0.025, 0.975])
            contrasts.append(
                dict(
                    panel=panel,
                    cap=cap,
                    contrast=col,
                    n=len(rows),
                    mean=x.mean(),
                    median=x.median(),
                    negative=(x < -1e-10).sum(),
                    unchanged=(abs(x) <= 1e-10).sum(),
                    positive=(x > 1e-10).sum(),
                    query_block_lo=limits[0],
                    query_block_hi=limits[1],
                )
            )
    for name, table in (
        ("distributions", pd.DataFrame(dist)),
        ("pairs", pairs),
        ("contrasts", pd.DataFrame(contrasts)),
    ):
        table.to_csv(directory / (name + ".csv"), index=False, float_format="%.9g")


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    summarize(p.parse_args().directory)
