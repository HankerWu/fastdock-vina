"""Plot matched released-pool accuracy at the two registered GBA1 search budgets."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent / "results"
fig, axes = plt.subplots(
    1, 2, figsize=(9, 3.65), sharex=True, sharey=True, layout="constrained"
)
labels = {"native": "FastDock Pro", "chpi": "Same poses, CH–π ranking"}
for ax, d, cap in zip(axes, (root, root / "cap24000"), ("12,000", "24,000")):
    with (d / "cells.csv").open() as f:
        rows = [r for r in csv.DictReader(f) if r["state"] == "protonated"]
    for arm, color, linestyle in (
        ("native", "#426578", "-"),
        ("chpi", "#b24322", "--"),
    ):
        x = sorted(float(r["top1"]) for r in rows if r["arm"] == arm)
        ax.step(
            [0, *x, 13.5],
            [0, *np.arange(1, len(x) + 1) / len(x), 1],
            where="post",
            color=color,
            linestyle=linestyle,
            linewidth=2,
            label=labels[arm],
        )
    ax.set_title(f"{cap} evaluations per trajectory", loc="left", fontsize=11)
    ax.set_xlabel("Top-1 heavy-atom RMSD (Å)")
    ax.grid(axis="y", alpha=0.2)
    ax.set(xlim=(0, 13.5), ylim=(0, 1.02))
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("Fraction of 27 conditions")
axes[1].legend(loc="lower right", frameon=False, fontsize=8)
fig.savefig(root / "top1_budget_ecdf.png", dpi=160)
