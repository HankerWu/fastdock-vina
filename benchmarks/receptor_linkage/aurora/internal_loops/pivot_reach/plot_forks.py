"""Show chronological local outcomes and actual released modes after each single injection."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap

ROOT = Path(__file__).resolve().parents[5]
HERE = Path(__file__).resolve().parent
data = pd.read_csv(HERE.parent / "rigid_capture/results/branches.csv")
columns = ["start_A", "after_hunt_A", "end_trajectory_min_A", "oracle_A", "top1_A"]
labels = ["Placed", "Hunt", "Trajectory\nend*", "Pool\nbest", "Rank 1"]
palette = ["#afdccc", "#dfede8", "#f2ede7", "#eed4bd", "#d8a889"]
bounds = [0, 2.000000001, 3.000000001, 5.000000001, 8.000000001, 15]
fig, axes = plt.subplots(1, 2, figsize=(9, 5.8), sharey=True)
for tr, ax in enumerate(axes, 1):
    subset = data[data.trajectory == tr].sort_values("step")
    matrix = subset[columns].to_numpy()
    ax.imshow(
        matrix,
        aspect="auto",
        cmap=ListedColormap(palette),
        norm=BoundaryNorm(bounds, len(palette)),
    )
    for (i, j), value in np.ndenumerate(matrix):
        ax.text(
            j,
            i,
            f"{value:.2f}",
            ha="center",
            va="center",
            fontsize=9,
            weight="bold" if value <= 2 else "normal",
            color="#162d36",
        )
    ax.set_xticks(range(len(labels)), labels)
    ax.set_yticks(range(len(subset)), subset.step)
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="both", length=0, pad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(f"Trajectory {tr}", loc="left", fontsize=12)
axes[0].set_ylabel("Step of single crystal-informed placement")
fig.suptitle(
    "AT9283  |  What survives after a near-crystal placement?",
    x=0.07,
    ha="left",
    fontsize=15,
)
fig.text(
    0.07,
    0.044,
    "Receptor-frame symmetry-corrected heavy-atom RMSD (Å).  Bold ≤2 Å.  *Best in the perturbed trajectory.",
    fontsize=9,
)
fig.text(
    0.07,
    0.014,
    "Each row restarts a separate docking run. Pool best and rank 1 include both trajectories.",
    fontsize=9,
)
fig.subplots_adjust(top=0.86, bottom=0.17, left=0.085, right=0.99, wspace=0.12)
fig.savefig(ROOT / "docs/figures/at9283_rigid_forks.png", dpi=155, bbox_inches="tight")
