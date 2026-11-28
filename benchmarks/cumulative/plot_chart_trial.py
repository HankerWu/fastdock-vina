"""ECDF for each individual fixed cross-docking cell and each registered BRD4 cell."""

import csv
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
COLORS = {"official": "#424e60", "stable": "#cf7641", "candidate": "#136c74"}
LABELS = {
    "official": "Official Vina / parent",
    "stable": "Starting Pro",
    "candidate": "Angle-limited trial",
}


def read(name):
    with (HERE.parent / name / "results/cells.csv").open() as f:
        return list(csv.DictReader(f))


def main():
    fixed = read("search_trial_panel")
    external = read("search_trial_nov28")
    fig, axes = plt.subplots(
        1, 2, figsize=(9.3, 3.9), sharey=True, layout="constrained"
    )
    for ax, rows, title in zip(
        axes,
        ([r for r in fixed if r["cohort"] == "reciprocal"], external),
        ("Fixed reciprocal panel · 20 cells", "BRD4 transfer · 8 cells"),
    ):
        arms = {}
        for arm in ("official", "parent", "stable", "candidate"):
            arms[arm] = sorted(
                float(r["first"])
                for r in rows
                if r["arm"] == arm and r["status"] == "ok"
            )
        assert (
            arms["official"] == arms["parent"]
            and len(arms["official"]) == len(rows) // 4
        )
        for arm, color in COLORS.items():
            x = arms[arm]
            ax.step(
                [0, *x, 17],
                [0, *np.arange(1, len(x) + 1) / len(x), 1],
                where="post",
                lw=2.2,
                color=color,
                label=LABELS[arm],
            )
        ax.grid(axis="y", color="#e2e7e8")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set(
            xlim=(0, 17),
            ylim=(0, 1.02),
            xticks=[0, 4, 8, 12, 16],
            xlabel="Top-1 heavy-atom RMSD (Å)",
        )
        ax.set_title(title, loc="left", fontsize=12)
    axes[0].set_ylabel("Fraction of cells at or below RMSD")
    axes[1].legend(loc="lower right", frameon=False, fontsize=9)
    fig.savefig(HERE / "chart_trial_ecdf.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
