"""Check whether unanimous seed/receptor agreement identifies trustworthy same-assay pairs."""

import argparse
import csv
import statistics
from collections import defaultdict
from decimal import Decimal
from pathlib import Path


def load(path):
    with path.open() as fp:
        return list(csv.DictReader(fp))


def save(path, rows):
    with path.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def compare(arms, output):
    records = []
    for (panel, arm, strong, weak), gaps in sorted(arms.items()):
        votes = list(gaps.values())
        expected = {"fxa": 4, "brd4": 12, "hcaii": 8, "t4l": 8}[panel]
        if len(votes) != expected:
            raise ValueError(f"Incomplete comparison for {panel} {strong}/{weak}")
        median = statistics.median(votes)
        correct = lambda x: 1 if x > 0 else 0 if x < 0 else 0.5
        records.append(
            dict(
                panel=panel,
                arm=arm,
                strong=strong,
                weak=weak,
                conditions=len(votes),
                votes_for_strong=sum(v > 0 for v in votes),
                votes_for_weak=sum(v < 0 for v in votes),
                median_gap=median,
                median_credit=correct(median),
                unanimous=int(all(v > 0 for v in votes) or all(v < 0 for v in votes)),
            )
        )
    save(output, records)
    for panel, arm in sorted({(r["panel"], r["arm"]) for r in records}):
        rows = [r for r in records if (r["panel"], r["arm"]) == (panel, arm)]
        stable = [r for r in rows if r["unanimous"]]
        print(
            panel,
            arm,
            "all",
            sum(r["median_credit"] for r in rows),
            len(rows),
            "unanimous",
            sum(r["median_credit"] for r in stable),
            len(stable),
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mapping_results", type=Path)
    p.add_argument("fxa_results", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    groups = defaultdict(dict)
    for panel in ("brd4", "t4l", "hcaii"):
        for r in load(a.mapping_results / panel / "pairs.csv"):
            k = (panel, r["arm"], r["strong"], r["weak"])
            condition = (r["receptor"], r["seed"])
            if condition in groups[k]:
                raise ValueError("Duplicate comparison")
            groups[k][condition] = round(Decimal(r["gap"]), 3)
    fxa = defaultdict(dict)
    for rec in ("1FJS", "2P16"):
        for r in load(a.fxa_results / f"fxa_{rec}_cells.csv"):
            fxa[(rec, r["seed"])][r["ligand"]] = r
    for condition, ligands in fxa.items():
        for s, v in ligands.items():
            for w, u in ligands.items():
                if Decimal(u["measurement"]) < 3 * Decimal(v["measurement"]):
                    continue
                for col in ("native", "no_divisor"):
                    groups["fxa", col, s, w][condition] = Decimal(u[col]) - Decimal(
                        v[col]
                    )
    compare(groups, a.output)
