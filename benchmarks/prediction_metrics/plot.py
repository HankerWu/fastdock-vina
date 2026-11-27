"""Plot two balanced held-out cross-docking strata at their original search budget."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "ligand_hydroxyl/transfer/results/cells.csv"


def plot(out):
    frame = pd.read_csv(SOURCE)
    frame = frame.loc[
        (frame.engine == "pro") & (frame.kind == "cross") & (frame.cap == 12000)
    ]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.1), sharey=True)
    styles = {
        "native": ("#304760", "-", "Native Pro"),
        "bonded": ("#ce5a30", "--", "Ligand OH integral"),
    }
    for ax, target, title in zip(
        axes, ("ERRg", "HSP90a"), ("ERRγ", "HSP90α"), strict=True
    ):
        pair = frame[frame.target.eq(target) & frame.arm.isin(styles)].pivot(
            index=["query", "receiver", "seed"], columns="arm", values="top1_A"
        )
        if (
            pair.isna().any().any()
            or not pair.groupby(level=["query", "receiver"]).size().eq(2).all()
        ):
            raise ValueError("Expected complete paired runs with two seeds per complex")
        for arm, (color, linestyle, label) in styles.items():
            values = np.sort(pair[arm].to_numpy())
            x = np.r_[0.0, values, 9.0]
            if x[-2] > x[-1]:
                raise ValueError("Displayed range would truncate the ECDF")
            y = np.r_[0.0, np.arange(1, len(values) + 1) / len(values), 1.0]
            ax.step(
                x,
                y,
                where="post",
                color=color,
                linewidth=2.5,
                linestyle=linestyle,
                label=label,
            )
        ax.axvline(2, color="#74797f", linewidth=0.9, linestyle=":")
        ax.text(2.09, 0.045, "2 Å", color="#676e76", fontsize=9, va="bottom")
        ax.set_title(
            f"{title}  ·  {len(pair) // 2} ligands, {len(pair)} runs",
            loc="left",
            fontsize=12,
        )
        ax.set_xlim(0, 9)
        ax.set_xticks([0, 2, 4, 6, 8])
        ax.set_ylim(-0.015, 1.04)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["bottom", "left"]].set_color("#849199")
        ax.grid(axis="y", color="#d9e0e3", linewidth=0.7)
        ax.set_axisbelow(True)
        ax.set_xlabel("Top-1 ligand RMSD threshold (Å)")
    axes[0].set_ylabel("Fraction of runs at or below threshold")
    axes[0].legend(loc="lower right", frameon=False, fontsize=9)
    fig.suptitle(
        "Cross-docking on the same candidate poses",
        x=0.08,
        y=0.98,
        ha="left",
        fontsize=16,
        fontweight="semibold",
    )
    fig.subplots_adjust(top=0.79, bottom=0.16, left=0.09, right=0.985, wspace=0.15)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160, bbox_inches="tight", pad_inches=0.13)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=HERE / "results/pose_ecdf.png")
    plot(parser.parse_args().out)
