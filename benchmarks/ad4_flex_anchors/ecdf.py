"""Generate the ECDF of actual first poses and their per-condition changes."""

import argparse
from pathlib import Path
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main(folder):
    with (folder / "cases.csv").open() as f:
        rows = list(csv.DictReader(f))
    keys = [(r["query"], r["receiver"], r["seed"]) for r in rows if r["arm"] == "old"]
    series = {}
    for arm in ["official", "old", "new"]:
        values = {
            (r["query"], r["receiver"], r["seed"]): float(r["top1"])
            for r in rows
            if r["arm"] == arm and r["status"] == "0" and r["top1"]
        }
        if len(keys) != len(set(keys)) or set(values) != set(keys):
            raise ValueError("Missing or duplicated structural comparisons")
        series[arm] = np.array([values[key] for key in keys])
    series["paired_delta"] = series["new"] - series["old"]
    fig, axes = plt.subplots(
        1, 2, figsize=(9.4, 3.6), layout="constrained", sharey=True
    )
    labels = {
        "official": "Official Vina 1.2.7 (AD4)",
        "old": "Pro before fix",
        "new": "Pro with anchor contacts",
        "paired_delta": "New minus old",
    }
    colors = {
        "official": "#697681",
        "old": "#bf6b37",
        "new": "#116b85",
        "paired_delta": "#116b85",
    }
    data = []
    for arm, values in series.items():
        x = np.sort(values)
        y = np.arange(1, len(x) + 1) / len(x)
        for value, count in zip(*np.unique(x, return_counts=True)):
            data.append(
                dict(
                    curve=arm,
                    distance_A=value,
                    count=count,
                    n=len(x),
                    cumulative=np.mean(x <= value),
                )
            )
        ax = axes[int(arm == "paired_delta")]
        lower = min(0.0, x[0] - 0.2)
        ax.step(
            np.r_[lower, x],
            np.r_[0.0, y],
            where="post",
            color=colors[arm],
            lw=2,
            label=labels[arm],
        )
    axes[0].legend(loc="lower right", frameon=False, fontsize=8)
    axes[0].set(
        xlabel="Top-1 heavy-atom RMSD (Å)",
        title="First-ranked pose",
        ylabel="Fraction of paired conditions",
    )
    axes[1].axvline(0, color="#667681", ls=":", lw=1)
    axes[1].set(xlabel="Change in Top-1 RMSD (Å)", title="Within-condition change")
    for ax in axes:
        ax.grid(axis="y", alpha=0.2)
        ax.set_ylim(0, 1.025)
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(folder / "top1_ecdf.png", dpi=150)
    plt.close(fig)
    with (folder / "ecdf.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, data[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    main(p.parse_args().folder)
