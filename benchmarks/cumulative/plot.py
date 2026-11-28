"""Plot every fixed-budget Top-1 observation in its own docking cohort."""

from pathlib import Path
import csv
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "geometry_comparison/upstream/four_versions/results/cells.csv"


def main():
    with DATA.open() as stream:
        rows = list(csv.DictReader(stream))
    fig, axes = plt.subplots(
        1, 2, figsize=(10, 4), sharey=True, constrained_layout=True
    )
    for ax, cohort, title in zip(
        axes, ("cognate", "reciprocal"), ("Cognate docking", "Reciprocal cross-docking")
    ):

        def measurements(arm):
            records = sorted(
                (r["query"], r["receiver"], r["seed"], float(r["first"]))
                for r in rows
                if r["arm"] == arm and r["cohort"] == cohort and r["status"] == "ok"
            )
            if len(records) != 20:
                raise ValueError("wrong cohort size")
            return records

        for arm, other, label, color in (
            ("official", "parent", "Official Vina and parent FastDock", "#276480"),
            ("candidate", "stable", "Initial and final Pro", "#AD4C32"),
        ):
            data = measurements(arm)
            if data != measurements(other):
                raise ValueError(
                    "versions differ; render separate ECDFs before updating the figure"
                )
            x = sorted(row[3] for row in data)
            ax.step(
                [0, *x, 17],
                [0, *np.arange(1, 21) / 20, 1],
                where="post",
                lw=2.25,
                color=color,
                label=label,
            )
        ax.set_title(title, loc="left", fontweight="bold", pad=10)
        ax.set(
            xlabel="Top-1 heavy-atom RMSD (Å)",
            xlim=(0, 17),
            ylim=(0, 1.02),
            xticks=range(0, 17, 4),
        )
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#d9dfe1", lw=0.7)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Fraction at or below the stated RMSD")
    axes[1].legend(loc="lower right", fontsize=8, frameon=True, framealpha=0.96)
    fig.savefig(HERE / "fixed_pose_ecdf.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
