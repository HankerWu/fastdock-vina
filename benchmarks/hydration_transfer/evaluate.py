"""Replay frozen water scoring and summarize paired predictions and oracle-start diagnostics."""

import argparse
import csv
import json
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from panel import HERE, PANELS, write_csv
from scipy.spatial import cKDTree

sys.path.insert(0, str(HERE.parent / "mobile_water"))
from mobile_water import Site
from replay import cross_sites, models, pdb_lines, xyz


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def evaluate(work, output):
    output.mkdir(parents=True, exist_ok=True)
    meta = json.loads((work / "metadata.json").read_text())
    geometry = {
        (
            r["query"],
            r["receiver"],
            r["cap"],
            r["seed"],
            r["engine"],
            int(r["model"]),
        ): r
        for r in rows(work / "geometry.csv")
    }
    costs, all_scores, cells, node_cost, site_rows = [], [], [], [], []
    summary = defaultdict(list)
    for entries in PANELS.values():
        receiver = next(iter(entries))
        t0 = time.perf_counter()
        pdb = pdb_lines(work / "raw" / (receiver + ".pdb"))
        protein = cKDTree(
            [xyz(l) for l in pdb if l.startswith("ATOM") and l[21] == "A"]
        )
        sites = cross_sites(pdb, protein, np.asarray(meta[receiver]["center"]), "A")
        select_s = time.perf_counter() - t0
        for label, center, b_iso in sites:
            site_rows.append(
                {
                    "receiver": receiver,
                    "label": label,
                    "x": center[0],
                    "y": center[1],
                    "z": center[2],
                    "B": b_iso,
                }
            )
        grids = {}
        for order in (11, 17):
            then = time.perf_counter()
            grids[order] = [
                Site.prepare(center, b, protein, order) for _, center, b in sites
            ]
            node_cost.append(
                {
                    "receiver": receiver,
                    "n_sites": len(sites),
                    "order": order,
                    "select_s": select_s,
                    "prepare_s": time.perf_counter() - then,
                }
            )
        for path in sorted((work / "runs").glob(f"*_{receiver}_*.pdbqt")):
            query, rec, cap, seed, engine = path.stem.split("_")
            context = {
                "target": meta[query]["target"],
                "kind": "cognate" if query == rec else "cross",
                "query": query,
                "receiver": rec,
                "cap": int(cap),
                "seed": int(seed),
                "engine": engine,
            }
            pool = []
            blocks = path.read_text().split("ENDMDL")
            for model, dry, tree in models(path):
                geom = geometry[query, rec, cap, seed, engine, model]
                if float(geom["dry"]) != dry:
                    raise ValueError(f"Dry score mismatch in {path}")
                base = dict(
                    **context,
                    model=model,
                    dry=dry,
                    rmsd=float(geom["rmsd"]),
                    native_total=float(
                        next(
                            l.split()[-1]
                            for l in blocks[model - 1].splitlines()
                            if l.startswith("REMARK INTER + INTRA:")
                        )
                    ),
                )
                if engine == "pro":
                    for order, grid in grids.items():
                        start = time.perf_counter()
                        terms = [site.evaluate(tree) for site in grid]
                        costs.append(
                            dict(
                                **context,
                                model=model,
                                order=order,
                                seconds=time.perf_counter() - start,
                            )
                        )
                        if order == 11:
                            base["point"] = dry + sum(t[0] for t in terms)
                        base["mobile" if order == 11 else "order17"] = dry + sum(
                            t[1] for t in terms
                        )
                pool.append(base)
                all_scores.append(base)
            for arm in (
                ("dry", "point", "mobile", "order17") if engine == "pro" else ("dry",)
            ):
                ranked = sorted(pool, key=lambda b: (b[arm], b["model"]))
                rr = [r["rmsd"] for r in ranked]
                cell = dict(
                    **context,
                    arm=arm,
                    count=len(rr),
                    first=ranked[0]["model"],
                    top1_rmsd=rr[0],
                    top1=int(rr[0] <= 2),
                    top5=int(min(rr[:5]) <= 2),
                    oracle=int(min(rr) <= 2),
                    oracle_rmsd=min(rr),
                )
                cells.append(cell)
                summary[context["target"], context["kind"], cap, engine, arm].append(
                    cell
                )
    write_csv(output / "sites.csv", site_rows)
    write_csv(output / "receivers.csv", node_cost)
    write_csv(output / "rescore_cost.csv", costs)
    write_csv(
        output / "poses.csv",
        [
            {**dict.fromkeys(["point", "mobile", "order17"], ""), **r}
            for r in all_scores
        ],
    )
    write_csv(output / "cells.csv", cells)
    agg = []
    for (target, kind, cap, engine, arm), group in sorted(summary.items()):
        count = {k: sum(i[k] for i in group) for k in ("top1", "top5", "oracle")}
        agg.append(
            dict(
                target=target,
                kind=kind,
                cap=cap,
                engine=engine,
                arm=arm,
                n=len(group),
                **count,
                conditional_top1=(
                    count["top1"] / count["oracle"] if count["oracle"] else ""
                ),
                median_top1_rmsd=statistics.median(r["top1_rmsd"] for r in group),
            )
        )
    write_csv(output / "summary.csv", agg)
    # Compare the local result with each cell under one common, unscaled energy zero.
    diagnostic = []
    for row in rows(work / "crystal_starts.csv"):
        if row["arm"] != "best_projected":
            continue
        for cell in cells:
            if (
                (cell["query"], cell["receiver"]) != (row["source"], row["receiver"])
                or cell["arm"] != "dry"
                or cell["engine"] != "pro"
            ):
                continue
            comparable = [
                r["native_total"]
                for r in all_scores
                if all(
                    r[k] == cell[k]
                    for k in ("query", "receiver", "cap", "seed", "engine")
                )
            ]
            oracle_energy = float(row["inter_plus_intra"])
            diagnostic.append(
                dict(
                    **cell,
                    projected_initial_A=row["initial_rmsd"],
                    projected_local_A=row["rmsd"],
                    local_native_total=oracle_energy,
                    best_pool_total=min(comparable),
                    local_minus_best_pool=oracle_energy - min(comparable),
                )
            )
    write_csv(output / "oracle_starts.csv", diagnostic)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    evaluate(args.work, args.output)
