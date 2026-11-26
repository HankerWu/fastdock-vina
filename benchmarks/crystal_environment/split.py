"""Report fixed-pool selections separately for cross-docking and cognate receptors."""

from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path

from evaluate import read, save


def stratify(results, t4_runs):
    rows = []
    for r in read(results / "cells.csv"):
        rows.append((r, r["method"], r["hit"], r["oracle"]))
    t4_scores = read(t4_runs / "modes.csv")
    blocks = defaultdict(list)
    for r in t4_scores:
        blocks[r["receptor"], r["ligand"], r["seed"], r["arm"]].append(r)
    for group in blocks.values():
        best = min(group, key=lambda r: int(r["model"]))
        hit = lambda r: float(r["rmsd"]) <= 2 and int(r["sign_ok"])
        r = dict(best, target="t4l", pool=best["arm"])
        rows.append((r, "native", hit(best), int(any(hit(s) for s in group))))
    for r in read(results / "density_cells.csv"):
        rows.append((r, r["measure"], r["joint"], None))
    for r in read(results / "shape_cells.csv"):
        rows.append(
            (
                r,
                "_".join(r[k] for k in ("kind", "metric", "selection")),
                r["joint"],
                None,
            )
        )
    strata = defaultdict(list)
    for r, method, hit, oracle in rows:
        stratum = "cognate" if r["receptor"] == r["ligand"] else "cross"
        for query in ("all", r["ligand"]):
            for label in ("all", stratum):
                strata[r["target"], r["pool"], method, query, label].append(
                    (int(hit), oracle)
                )
    save(
        results / "splits.csv",
        [
            dict(
                target=t,
                pool=p,
                method=m,
                query=q,
                receptor_relation=s,
                n=len(values),
                hits=sum(v[0] for v in values),
                oracle=sum(int(v[1]) for v in values)
                if values[0][1] is not None
                else "",
            )
            for (t, p, m, q, s), values in sorted(strata.items())
        ],
    )


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("results", type=Path)
    p.add_argument("t4_runs", type=Path)
    a = p.parse_args()
    stratify(a.results, a.t4_runs)
