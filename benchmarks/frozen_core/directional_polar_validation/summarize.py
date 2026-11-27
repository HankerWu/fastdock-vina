"""Summarize all individual seed/budget cells and audit fixed-pool rank invariance."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path
import re


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def save(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(work):
    cells = read(work / "cells.csv")
    groups = defaultdict(list)
    natives = {}
    for c in cells:
        groups[c["pdb"], c["engine"], c["cap"], c["arm"]].append(c)
        if c["arm"] == "native":
            natives[c["pdb"], c["engine"], c["cap"], c["seed"]] = c
    rows = []
    for (pdb, engine, cap, arm), group in groups.items():
        bases = [natives[pdb, engine, cap, c["seed"]] for c in group]
        assert all(c["rmsd_oracle"] == b["rmsd_oracle"] for c,b in zip(group,bases))
        oracle = sum(int(c["oracle"]) for c in group)
        hits = sum(int(c["success1"]) for c in group)
        rows.append(dict(pdb=pdb, engine=engine, cap=cap, arm=arm, runs=len(group), top1=hits,
                         top5=sum(int(c["success5"]) for c in group), oracle=oracle,
                         conditional_top1=hits / oracle if oracle else "",
                         top1_mode_changes=sum(c["top1_mode"] != b["top1_mode"] for c,b in zip(group,bases))))
    save(work / "summary.csv", rows)

    counts = []
    for raw in sorted((work / "runs").glob("*_pro.log")):
        trajectories = [s for s in raw.read_text().splitlines() if s.startswith("Search trajectory:")]
        assert len(trajectories) == 4
        totals = {k: sum(int(re.search(r"\b" + k + r"=(\d+)", s)[1]) for s in trajectories)
                  for k in ("steps", "hunt_calls", "hunt_evals", "second_calls", "second_evals")}
        pdb, cap, seed, _ = raw.stem.split("_")
        explicit = int(re.search(r"Search counters: .*\bexplicit_evals=(\d+)", raw.read_text())[1])
        counts.append(dict(pdb=pdb, cap=cap, seed=seed, **totals, explicit_evals=explicit,
                           total_observed=totals["hunt_evals"]+totals["second_evals"]+explicit))
    save(work / "evaluation_counts.csv", counts)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    summarize(p.parse_args().work)
