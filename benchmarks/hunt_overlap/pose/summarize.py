"""Descriptive paired endpoints; reuse of a structure or seed does not create new targets."""

import argparse
import csv
import re
from pathlib import Path
from statistics import median


def load(path):
    with path.open() as file:
        return list(csv.DictReader(file))


def write(path, rows):
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    work = ap.parse_args().work
    rows = load(work / "cells.csv")
    groups = [("overall", "all", rows)]
    for name in ("target", "cap", "kind", "query"):
        groups.extend(
            (name, key, [r for r in rows if r[name] == key])
            for key in sorted({r[name] for r in rows})
        )
    result = []
    for field, level, group in groups:
        for arm in (
            "upstream",
            "parent",
            *sorted({r["arm"] for r in rows} - {"upstream", "parent"}),
        ):
            g = [r for r in group if r["arm"] == arm]
            hits = {
                key: sum(float(r[key]) <= 2 for r in g)
                for key in ("top1_A", "best5_A", "oracle_A")
            }
            energy_grad = sum(
                int(r[k])
                for r in g
                for k in ("hunt_evals", "second_evals", "explicit_evals")
            )
            calls = sum(
                int(r[k])
                for r in g
                for k in ("hunt_calls", "second_calls", "explicit_calls")
            )
            result.append(
                {
                    "group": field,
                    "level": level,
                    "arm": arm,
                    "n": len(g),
                    "top1": hits["top1_A"],
                    "top5": hits["best5_A"],
                    "oracle": hits["oracle_A"],
                    "conditional_top1": hits["top1_A"] / hits["oracle_A"]
                    if hits["oracle_A"]
                    else "",
                    "median_top1_A": median(float(r["top1_A"]) for r in g),
                    "median_oracle_A": median(float(r["oracle_A"]) for r in g),
                    "top1_with_clash": sum(int(r["top1_clash_pairs"]) > 0 for r in g),
                    "total_energy_gradient": energy_grad,
                    "total_optimizations": calls,
                    "total_seconds": sum(float(r["seconds"]) for r in g),
                }
            )
    write(work / "summary.csv", result)
    controls = work / "crystal_starts.csv"
    if controls.exists():
        guides = {
            (r["source"], r["receiver"]): r
            for r in load(controls)
            if r["arm"] == "best_projected"
        }
        comparisons = []
        for row in rows:
            guide = guides.get((row["query"], row["receiver"]))
            if row["arm"] != "parent" or guide is None:
                continue
            text = (work / "runs" / (row["key"] + ".pdbqt")).read_text()
            values = [
                float(x)
                for x in re.findall(
                    r"^REMARK INTER \+ INTRA:\s*(\S+)", text, re.MULTILINE
                )
            ]
            if len(values) != int(row["modes"]):
                raise ValueError("Missing energy component in " + row["key"])
            ideal = float(guide["inter_plus_intra"])
            comparisons.append(
                {
                    "key": row["key"],
                    "target": row["target"],
                    "pool_oracle_A": row["oracle_A"],
                    "best_projected_after_A": guide["rmsd"],
                    "lowest_pool_inter_plus_intra": min(values),
                    "guided_inter_plus_intra": ideal,
                    "guided_minus_best_pool": ideal - min(values),
                    # Inputs are printed at 0.001 kcal/mol; count possible ties conservatively.
                    "returned_modes": len(values),
                    "guided_rank_best": 1 + sum(e < ideal - 0.002 for e in values),
                    "guided_rank_worst": 1 + sum(e <= ideal + 0.002 for e in values),
                }
            )
        if comparisons:
            write(work / "guided_scores.csv", comparisons)
    print(*result[:4], sep="\n")


if __name__ == "__main__":
    main()
