"""Plot paired empirical distributions of released Aurora A poses."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parent / "results"
STYLES = {
    "official": ("Official Vina", "#38526b", "-"),
    "pro": ("Pro", "#c45032", "--"),
    "double": ("Pro, double budget", "#247568", ":"),
}


def plot():
    with (ROOT / "geometry_cells.csv").open() as source:
        cases = [r for r in csv.DictReader(source) if r["primary"] == "True"]
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.9), sharey=True)
    for ax, metric, title in zip(
        axes, ("top1", "oracle"), ("First-ranked pose", "Best released pose")
    ):
        for method, (label, color, line) in STYLES.items():
            data = sorted(float(r[metric]) for r in cases if r["method"] == method)
            if len(data) != 16 or data[-1] > 14:
                raise ValueError("Unexpected paired coverage or plot range")
            ax.step(
                [0, *data, 14],
                [0, *np.arange(1, 17) / 16, 1],
                where="post",
                label=label,
                color=color,
                linestyle=line,
                linewidth=2.2,
            )
        ax.set(xlim=(0, 14), ylim=(-0.025, 1.04), xlabel="Heavy-atom RMSD (Å)")
        ax.set_title(title, loc="left", fontsize=12, fontweight="semibold")
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.set_yticks([0, 0.25, 0.50, 0.75, 1])
        ax.grid(axis="y", linewidth=0.6, color="#d9dfe3")
        ax.spines[["right", "top"]].set_visible(False)
    axes[0].set_ylabel("Fraction of paired conditions")
    axes[0].legend(frameon=False, loc="upper left", fontsize=8.5)
    fig.suptitle(
        "Aurora A cross-docking  ·  two ligands, two receptors, four seeds",
        x=0.06,
        ha="left",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(ROOT / "pose_ecdf.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    plot()
