"""Combine verified classifications and plot the empirical receptor-frame RMSD CDFs."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from study import KEY


def main(results):
    keys = ["order", *KEY]
    pools = pd.read_csv(results / "cells.csv.gz", float_precision="round_trip")
    game = pd.read_csv(results / "game_pools.csv", float_precision="round_trip")
    combined = pools.merge(game[keys + ["lower", "upper"]], on=keys, how="left", validate="one_to_one")
    if not np.array_equal(combined.lower.isna(), combined.closed):
        raise ValueError("Game checks missing or unexpected")
    low = combined.lower.fillna(combined.best_possible_lower)
    high = combined.upper.fillna(combined.witness_upper)
    if np.max(np.abs(low-high)) > 1e-10:
        raise ArithmeticError("Open bounds cannot be presented as optimal values")
    combined["best_attainable"] = high
    combined["irreducible_gap"] = high - combined.oracle
    combined["dry_to_best"] = combined.dry_top - high
    combined.to_csv(results / "classified_pools.csv.gz", index=False,
                    compression={"method": "gzip", "mtime": 0})
    primary = combined[combined.order == 11]
    fig, ax = plt.subplots(figsize=(8, 4.8), layout="constrained")
    for column, label, color, style in (
            ("oracle","Best member of pool","#456c33","-"),
            ("best_attainable","Best selectable by any water law (uses labels)","#21748a","-"),
            ("dry_top","Dry top pose","#6450a3","--"),
            ("water_top","Fixed water law, p = 0.75","#ad4f29","-.")):
        x = np.sort(primary[column])
        ax.step(np.r_[0, x], np.r_[0, np.arange(1, len(x)+1)/len(x)], where="post",
                label=label, color=color, linestyle=style, linewidth=1.8)
    ax.set(xlim=(0, 15), ylim=(0, 1.015), xlabel="Receptor-frame heavy-atom RMSD (Å)",
           ylabel="Cumulative fraction of 156 pools")
    ax.spines[["right", "top"]].set_visible(False)
    ax.grid(axis="y", color="0.91")
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    fig.savefig(results / "distribution.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results",type=Path)
    main(parser.parse_args().results)
