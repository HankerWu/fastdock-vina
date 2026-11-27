"""Plot the shared-prior sensitivity of the conditional pose comparison."""
import csv
import math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def counts(rows, x):
    wins = np.zeros_like(x, dtype=int)
    for r in rows:
        a, b = int(r["sign_2J4Z"]), int(r["sign_2X6E"])
        c = r["crossover_log_odds_2X6E"]
        if c:
            z = float(c) / math.log(10)
            wins += (x < z) if a > 0 else (x > z)
        elif max(a, b) > 0:
            wins += 1
    return wins


def main():
    with (HERE / "results/pairs.csv").open() as f:
        rows = [r for r in csv.DictReader(f) if r["nodes"] == "65536"]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8), sharey=True, constrained_layout=True)
    for ax, frame, title in zip(axes, ("pocket", "whole"),
                                ("Pocket Cα registration", "Whole-chain Cα registration")):
        for field, color, label in (("center", "#747b85", "Fixed center"),
                                    ("window", "#006e85", "Local rigid-motion integral")):
            sub = [r for r in rows if r["registration"] == frame and r["model"] == field]
            breaks = sorted({float(r["crossover_log_odds_2X6E"]) / math.log(10)
                             for r in sub if r["crossover_log_odds_2X6E"]})
            edges = np.array([-10., *[b for b in breaks if -10 < b < 10], 10.])
            values = counts(sub, (edges[1:]+edges[:-1])/2)
            ax.stairs(values, edges, baseline=None, lw=2, color=color, label=label)
            ax.plot(0, counts(sub, np.array([0.]))[0], "o", color=color, ms=4)
        ax.axvline(0, ls=":", lw=1, color="#666666", zorder=0)
        ax.set(title=title, xlabel="log₁₀ prior odds of 2X6E vs 2J4Z", xlim=(-10, 10),
               ylim=(-.4,16.4), yticks=[0,2,4,6,8,10,12,14,16])
        ax.grid(axis="y", alpha=.18)
        ax.spines[["top","right"]].set_visible(False)
    axes[0].set_ylabel("Crystal-guided witness preferred (of 16)")
    axes[1].legend(loc="upper center", frameon=False, fontsize=9)
    fig.savefig(HERE / "results/prior_response.png", dpi=160)


if __name__ == "__main__":
    main()
