"""Summarize paired success, candidate-pool ceilings, RMSD and cost."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from statistics import median


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def seed_cohort(seed):
    value = int(seed)
    if value in (51913, 51941, 51973):
        return "original"
    if value in (52121, 52147, 52163):
        return "fresh"
    if value in (62347, 62383):
        return "deep"
    raise ValueError(f"Unregistered seed {value}")


def summary(work):
    by_run = defaultdict(list)
    for p in sorted(work.glob("modes-*.csv")):
        for r in csv.DictReader(p.open()):
            by_run[
                tuple(
                    r[k]
                    for k in (
                        "cap",
                        "target",
                        "ligand",
                        "receiver",
                        "kind",
                        "seed",
                        "arm",
                    )
                )
            ].append(r)
    cells = []
    for key, rows in sorted(by_run.items()):
        first, *rest = rows
        cells.append(
            dict(
                zip(
                    ("cap", "target", "ligand", "receiver", "kind", "seed", "arm"), key
                ),
                top1=float(first["rmsd"]),
                oracle=min(float(r["rmsd"]) for r in rows),
                modes=len(rows),
                evals=int(first["hunt"]) + int(first["second"]),
            )
        )
    write(work / "cells.csv", cells)
    report = []
    for cap in sorted({r["cap"] for r in cells}, key=int):
        for candidate in sorted(
            {r["arm"] for r in cells if r["cap"] == cap} - {"native"}
        ):
            for cohort in sorted(
                {
                    seed_cohort(r["seed"])
                    for r in cells
                    if r["cap"] == cap and r["arm"] == candidate
                }
            ):
                for name in [
                    "overall",
                    "cognate",
                    "cross",
                    *sorted({r["target"] for r in cells}),
                ]:
                    g = [
                        r
                        for r in cells
                        if r["cap"] == cap
                        and (seed_cohort(r["seed"]) == cohort)
                        and name in ("overall", r["target"], r["kind"])
                    ]
                    ar = [
                        [x for x in g if x["arm"] == arm]
                        for arm in ["native", candidate]
                    ]
                    if not ar[0] and not ar[1]:
                        continue
                    assert len(ar[0]) == len(ar[1])
                    assert [(q["ligand"], q["receiver"], q["seed"]) for q in ar[0]] == [
                        (q["ligand"], q["receiver"], q["seed"]) for q in ar[1]
                    ]
                    d = dict(
                        cap=cap,
                        cohort=cohort,
                        group=name,
                        candidate=candidate,
                        n=len(ar[0]),
                    )
                    for key in ["top1", "oracle"]:
                        for a, arm in zip(ar, ("native", candidate)):
                            d[
                                ("native" if arm == "native" else "candidate")
                                + "_"
                                + key
                            ] = sum(x[key] <= 2 for x in a)
                            d[
                                ("native" if arm == "native" else "candidate")
                                + "_"
                                + key
                                + "_median"
                            ] = median(x[key] for x in a)
                        d[key + "_gain"] = sum(
                            b[key] <= 2 < a[key] for a, b in zip(*ar)
                        )
                        d[key + "_loss"] = sum(
                            a[key] <= 2 < b[key] for a, b in zip(*ar)
                        )
                    for a, arm in zip(ar, ("native", candidate)):
                        d[
                            ("native" if arm == "native" else "candidate")
                            + "_mean_evals"
                        ] = sum(x["evals"] for x in a) / len(a)
                    report.append(d)
    write(work / "summary.csv", report)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    a = p.parse_args()
    summary(a.work)
