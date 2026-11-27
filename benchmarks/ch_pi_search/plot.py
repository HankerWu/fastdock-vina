"""Plot top-pose empirical cumulative distributions on the two iminosugar panels."""

from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot(data, path):
    colors = ["#545b61", "#267968", "#c06036", "#424ba5"]
    configs = [
        ("native", "native", "Baseline"),
        ("field", "native", "Search only"),
        ("native", "field", "Rerank only"),
        ("field", "field", "Search + rerank"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.4), sharex=True, sharey=True)
    for i, panel in enumerate(["GBA1", "GAA"]):
        for j, cap in enumerate([12000, 24000]):
            ax = axes[i, j]
            crop = data[(data.panel == panel) & (data.cap == cap)]
            for (s, r, label), color in zip(configs, colors):
                values = np.sort(crop[(crop.search == s) & (crop.rerank == r)].top1)
                x = np.r_[0, values, 14]
                y = np.r_[0, np.arange(1, len(values) + 1) / len(values), 1]
                ax.step(x, y, where="post", label=label, c=color, lw=1.8)
            ax.set_title(
                f"{panel}  ·  {cap:,} evaluations per trajectory",
                loc="left",
                fontsize=10,
            )
            ax.set_xlim(0, 14)
            ax.set_ylim(-0.025, 1.04)
            ax.grid(alpha=0.18)
            ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("Fraction of conditions")
    for ax in axes[-1]:
        ax.set_xlabel("Top-1 heavy-atom RMSD (Å)")
    fig.legend(
        *axes[0, 0].get_legend_handles_labels(),
        loc="upper center",
        ncol=4,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("history", type=Path)
    p.add_argument("gaa", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    plot(pd.concat([pd.read_csv(a.history), pd.read_csv(a.gaa)]), a.out)
