"""Paired continuous errors by ligand and seed batch; one observation per run."""

import argparse
import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
FIELDS = ("top1_A", "best5_A", "oracle_A")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=Path, default=ROOT / "results")
    parser.add_argument("--repeat", type=Path, default=ROOT / "replication")
    parser.add_argument("--output", type=Path, default=ROOT / "summary.csv")
    args = parser.parse_args()
    rows = []
    for batch, folder in (("results", args.first), ("replication", args.repeat)):
        with (folder / "cells.csv").open() as stream:
            cells = list(csv.DictReader(stream))
        for q in ("all", "4EY6", "4EY7"):
            groups = {}
            for state in ("neutral", "cation"):
                g = sorted(
                    (
                        r
                        for r in cells
                        if r["state"] == state and q in ("all", r["query"])
                    ),
                    key=lambda r: (r["query"], r["receiver"], r["seed"]),
                )
                groups[state] = g

            def keys(group):
                return [(r["query"], r["receiver"], r["seed"]) for r in group]

            assert keys(groups["neutral"]) == keys(groups["cation"])
            for metric in FIELDS:
                v = {
                    k: np.array([float(r[metric]) for r in g])
                    for k, g in groups.items()
                }
                v["paired_delta"] = v["cation"] - v["neutral"]
                for key, x in v.items():
                    quant = np.quantile(
                        x, [0.25, 0.5, 0.75, 0.9, 0.95], method="averaged_inverted_cdf"
                    )
                    rows.append(
                        dict(
                            batch=batch,
                            query=q,
                            metric=metric,
                            state=key,
                            n=len(x),
                            mean=float(np.mean(x)),
                            q25=quant[0],
                            median=quant[1],
                            q75=quant[2],
                            p90=quant[3],
                            p95=quant[4],
                            n_lower=int((x < 0).sum()) if key == "paired_delta" else "",
                        )
                    )
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
