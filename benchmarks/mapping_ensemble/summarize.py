"""Summarize paired pose hits and independent chemical pairs by seed split."""

import argparse
from collections import defaultdict
from pathlib import Path
from statistics import median

from pair_stability import load, save


def summaries(folder, panel):
    cells = [r for r in load(folder / "cells.csv") if r["method"] == "released"]
    pairs = load(folder / "pairs.csv")
    for split, seeds in (
        ("first", {"314159", "271828"}),
        ("fresh", {"161803", "141421"}),
        ("all", {"314159", "271828", "161803", "141421"}),
    ):
        for arm in ("native", "single", "multi"):
            chosen = [r for r in cells if r["arm"] == arm and r["seed"] in seeds]
            comparisons = [r for r in pairs if r["arm"] == arm and r["seed"] in seeds]
            gaps = defaultdict(list)
            for r in comparisons:
                gaps[r["strong"], r["weak"]].append(round(float(r["gap"]) * 1000))
            medians = [median(x) for x in gaps.values()]
            yield dict(
                panel=panel,
                split=split,
                arm=arm,
                cells=len(chosen),
                top=sum(int(r["hit"]) for r in chosen),
                signed=sum(int(r["strict"]) for r in chosen),
                oracle=sum(float(r["best_rmsd"]) <= 2 for r in chosen),
                median_rmsd=round(median(float(r["rmsd"]) for r in chosen), 5),
                repeated_credit=sum(float(r["credit"]) for r in comparisons),
                repeated_total=len(comparisons),
                unique_credit=sum(
                    1 if x > 0 else 0.5 if x == 0 else 0 for x in medians
                ),
                unique_total=len(medians),
            )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results", type=Path)
    a = p.parse_args()
    panels = ["t4l", "brd4", "hcaii"]
    if (a.results / "hcaii_declared").is_dir():
        panels.append("hcaii_declared")
    save(
        a.results / "summary.csv",
        [r for panel in panels for r in summaries(a.results / panel, panel)],
    )
