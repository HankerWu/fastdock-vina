"""Aggregate repeated seeds by cohort/budget and verify matched first-trial coordinates."""

import argparse
import csv
from pathlib import Path

from run_panel import write


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("panel", type=Path)
    p.add_argument("first", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    rows = [
        r for f in sorted(a.panel.glob("*_cells.csv")) for r in csv.DictReader(f.open())
    ]
    runs = []
    for split, cap, arm in sorted({(r["split"], r["cap"], r["arm"]) for r in rows}):
        rr = [r for r in rows if (r["split"], r["cap"], r["arm"]) == (split, cap, arm)]
        n = sum(int(r["top1"]) for r in rr)
        oracle = sum(int(r["oracle"]) for r in rr)
        runs.append(
            dict(
                split=split,
                cap=cap,
                arm=arm,
                cells=len(rr),
                top1=n,
                top5=sum(int(r["top5"]) for r in rr),
                oracle=oracle,
                conditional_top1=n / oracle if oracle else "",
                wall_s=sum(float(r["wall_s"]) for r in rr),
                cpu_s=sum(float(r["cpu_s"]) for r in rr),
                evals=sum(int(r["evals"]) for r in rr)
                if all(r["evals"] for r in rr)
                else "",
                peak_rss_KiB=max(int(r["peak_rss_KiB"]) for r in rr),
            )
        )
    write(a.out / "summary.csv", runs)
    trials = []
    for f in sorted(a.first.glob("*-parent.csv")):
        data = {
            arm: list(csv.DictReader(f.with_name(f.name.replace("parent", arm)).open()))
            for arm in ("parent", "wall", "joint")
        }
        for arm in ("wall", "joint"):
            assert [
                [r[k] for k in ("atom", "type", "charge", "x", "y", "z")]
                for r in data["parent"]
            ] == [
                [r[k] for k in ("atom", "type", "charge", "x", "y", "z")]
                for r in data[arm]
            ]
        e = {arm: float(d[0]["map_total"]) for arm, d in data.items()}

        def max_diff(a, b):
            return max(
                abs(float(x[key]) - float(y[key]))
                for x, y in zip(data[a], data[b])
                for key in ("gx", "gy", "gz")
            )

        pdb, seed, _ = f.stem.split("-")
        trials.append(
            dict(
                pdb=pdb,
                seed=seed,
                atoms=len(data["parent"]),
                identical_input=True,
                **e,
                wall_minus_parent=e["wall"] - e["parent"],
                joint_minus_wall=e["joint"] - e["wall"],
                max_abs_dgradient_wall=max_diff("parent", "wall"),
                max_abs_dgradient_charge=max_diff("wall", "joint"),
            )
        )
    write(a.out / "first_trials.csv", trials)


if __name__ == "__main__":
    main()
