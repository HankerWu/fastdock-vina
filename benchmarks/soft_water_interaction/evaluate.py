"""Score both searches on their exact released coordinates and compare the four arms."""

import argparse
import itertools
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
for sibling in ("mobile_water", "paired_crossdocking", "torsion_proposals"):
    sys.path.insert(0, str(HERE.parent / sibling))
from mobile_water import Site
from paired import aligned_receivers
from replay import cross_sites, models, pdb_lines, records, save, xyz
from run import alignment
from transfer import BUDGETS, SEEDS

KEY = ["panel", "target", "cap", "query", "receiver", "seed"]
ARMS = ["native_dry", "soft_dry", "native_water", "soft_water"]


def receptor_sites(work, receiver, center):
    lines = pdb_lines(work / "raw" / f"{receiver}.pdb")
    protein = cKDTree(
        [xyz(line) for line in lines if line.startswith("ATOM") and line[21] == "A"]
    )
    return [
        Site.prepare(pos, b, protein)
        for _, pos, b in cross_sites(lines, protein, np.array(center))
    ]


def historical_index():
    result = {}
    for cap in BUDGETS:
        table = HERE.parent / "paired_crossdocking" / "results" / f"poses_{cap}.csv"
        for row in records(table):
            if row["arm"] in ("native", "soft") and row["kind"] == "cross":
                result[
                    cap,
                    row["ligand"],
                    row["receiver"],
                    int(row["seed"]),
                    row["arm"],
                    int(row["rank"]),
                ] = float(row["rmsd"])
    old_water = {}
    for row in records(HERE.parent / "mobile_water/results/poses.csv"):
        if row["target"] != "thrombin":
            old_water[
                int(row["budget"]),
                row["ligand"],
                row["receiver"],
                int(row["seed"]),
                "native",
                int(row["model"]),
            ] = float(row["mobile"])
    return result, old_water


def collect(work, panel, history, old_water):
    meta = json.loads((work / "metadata.json").read_text())
    then = time.perf_counter()
    sites = {rec: receptor_sites(work, rec, m["center"]) for rec, m in meta.items()}
    site_seconds = time.perf_counter() - then
    rows, costs, seen = [], [], set()
    scoring_seconds = 0.0
    for q, r in itertools.permutations(sorted(meta), 2):
        if meta[q]["target"] != meta[r]["target"]:
            continue
        a, b = meta[q], meta[r]
        pair = aligned_receivers(a, b) if panel == "transfer" else (a, b)
        rot, trans, _, _ = alignment(*pair)
        refs = np.array(a["mapped"]) @ rot + trans
        for cap in BUDGETS:
            seeds = (
                SEEDS
                if panel == "transfer"
                else sorted({k[3] for k in history if k[:3] == (cap, q, r)})
            )
            for seed, search in itertools.product(seeds, ("native", "soft")):
                keys = dict(
                    panel=panel,
                    target=a["target"],
                    cap=cap,
                    query=q,
                    receiver=r,
                    seed=seed,
                    search=search,
                )
                path = work / f"runs_{cap}" / f"{q}_{r}_{seed}_{search}.pdbqt"
                payload = path.read_text().split("ENDMDL")
                if payload[-1].strip() or len(payload) < 2:
                    raise ValueError(f"Invalid released pool: {path}")
                for block, (model, dry, ligand) in zip(
                    payload[:-1], models(path), strict=True
                ):
                    then = time.perf_counter()
                    water = sum(site.evaluate(ligand)[1] for site in sites[r])
                    scoring_seconds += time.perf_counter() - then
                    atoms = {
                        int(line[6:11]): xyz(line)
                        for line in block.splitlines()
                        if line.startswith(("ATOM  ", "HETATM"))
                    }
                    coords = np.array([atoms[i] for i in a["serial"]])
                    rmsd = float(
                        np.sqrt(np.square(coords - refs).sum(-1).mean(-1)).min()
                    )
                    key = (cap, q, r, seed, search, model)
                    if panel == "historical":
                        seen.add(key)
                        if abs(history[key] - rmsd) > 1e-8:
                            raise ValueError(f"RMSD disagrees with original: {key}")
                        if (
                            search == "native"
                            and abs(old_water[key] - dry - water) > 1e-8
                        ):
                            raise ValueError(
                                f"Water score disagrees with original: {key}"
                            )
                    rows.append(
                        dict(
                            **keys,
                            model=model,
                            n_sites=len(sites[r]),
                            rmsd=rmsd,
                            dry=dry,
                            water_penalty=water,
                            hydrated=dry + water,
                        )
                    )
                log = path.with_suffix(".log").read_text()
                traj = re.findall(r"^Search trajectory:.*$", log, re.MULTILINE)
                if len(traj) != 4 or ("Soft receptor search:" in log) != (
                    search == "soft"
                ):
                    raise ValueError(f"Wrong trajectory metadata: {path}")
                sums = {
                    field: sum(int(re.search(rf"\b{field}=(\d+)", t)[1]) for t in traj)
                    for field in ("hunt_evals", "second_evals")
                }
                explicit = int(re.search(r"\bexplicit_evals=(\d+)", log)[1])
                costs.append(dict(**keys, **sums, explicit_evals=explicit))
    if panel == "historical" and (
        seen != history.keys() or not old_water.keys() <= seen
    ):
        raise ValueError("Historical cases missing")
    return (
        rows,
        costs,
        dict(
            panel=panel,
            prep_sites_s=site_seconds,
            evaluate_water_s=scoring_seconds,
            poses=len(rows),
        ),
    )


def describe(values):
    array = np.array(values)
    qs = np.quantile(array, [0.25, 0.50, 0.75, 0.9, 0.95])
    return dict(
        mean=array.mean(), q25=qs[0], median=qs[1], q75=qs[2], p90=qs[3], p95=qs[4]
    )


def summarize(raw, out):
    cells = []
    for ks, frame in raw.groupby(KEY + ["search"]):
        keys = dict(zip(KEY + ["search"], ks))
        for ranker, column in (("dry", "dry"), ("water", "hydrated")):
            ordered = frame.sort_values([column, "model"])
            top = ordered.iloc[0]
            if ranker == "dry" and top.model != 1:
                raise ValueError("Dry ranking fails released top-1")
            cells.append(
                dict(
                    **keys,
                    ranker=ranker,
                    arm=f"{keys['search']}_{ranker}",
                    n=len(frame),
                    n_sites=int(top.n_sites),
                    model=int(top.model),
                    rmsd=top.rmsd,
                    best5=ordered.rmsd.iloc[:5].min(),
                    oracle=ordered.rmsd.min(),
                    regret=top.rmsd - ordered.rmsd.min(),
                )
            )
    df = pd.DataFrame(cells)
    df.to_csv(out / "cells.csv", index=False)
    pairs = df.pivot(index=KEY, columns="arm", values="rmsd").reset_index()
    pairs.columns.name = None
    pairs["vs_baseline"] = pairs.soft_water - pairs.native_dry
    pairs["vs_soft"] = pairs.soft_water - pairs.soft_dry
    pairs["vs_water"] = pairs.soft_water - pairs.native_water
    pairs["interaction"] = pairs.vs_soft - (pairs.native_water - pairs.native_dry)
    pairs.to_csv(out / "pairs.csv", index=False)
    dist, deltas, intervals = [], [], []
    for panel, part in df.groupby("panel"):
        pairs_part = pairs.loc[pairs.panel == panel]
        groups = [("all", "all", part, pairs_part)]
        groups += [
            (t, "all", g, pairs_part.loc[pairs_part.target == t])
            for t, g in part.groupby("target")
        ]
        groups += [
            ("all", cap, g, pairs_part.loc[pairs_part.cap == cap])
            for cap, g in part.groupby("cap")
        ]
        for target, cap, group, paired in groups:
            k = dict(panel=panel, target=target, cap=cap)
            for arm, frame in group.groupby("arm"):
                for metric in ("rmsd", "best5", "oracle", "regret"):
                    dist.append(
                        dict(
                            **k,
                            arm=arm,
                            metric=metric,
                            n=len(frame),
                            **describe(frame[metric]),
                        )
                    )
            for metric in ("vs_baseline", "vs_soft", "vs_water", "interaction"):
                v = paired[metric]
                deltas.append(
                    dict(
                        **k,
                        metric=metric,
                        n=len(v),
                        better=(v < -1e-9).sum(),
                        worse=(v > 1e-9).sum(),
                        tied=(abs(v) <= 1e-9).sum(),
                        **describe(v),
                    )
                )
                if target == "all":
                    m = paired.groupby("target")[metric].mean()
                    draws = [
                        np.mean(d)
                        for d in itertools.product(m.to_numpy(), repeat=len(m))
                    ]
                    lo, hi = np.quantile(draws, [0.025, 0.975])
                    intervals.append(
                        dict(
                            **k,
                            metric=metric,
                            targets=len(m),
                            permutations=len(draws),
                            macro=m.mean(),
                            lo=lo,
                            hi=hi,
                        )
                    )
    for name, table in (
        ("distributions", dist),
        ("contrasts", deltas),
        ("target_bootstrap", intervals),
    ):
        save(out / f"{name}.csv", table)


def main(historical, transfer, out):
    out.mkdir(exist_ok=True, parents=True)
    history, old_water = historical_index()
    raw, timings, counters = [], [], []
    for work, panel in ((historical, "historical"), (transfer, "transfer")):
        rows, costs, timing = collect(work, panel, history, old_water)
        raw.extend(rows)
        counters.extend(costs)
        timings.append(timing)
    save(out / "poses.csv", raw)
    save(out / "counters.csv", counters)
    save(out / "scoring_cost.csv", timings)
    summarize(pd.DataFrame(raw), out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("historical", type=Path)
    p.add_argument("transfer", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    main(a.historical, a.transfer, a.out)
