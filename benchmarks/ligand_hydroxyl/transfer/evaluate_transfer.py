"""Score identical Pro candidates, measure endpoints and audit paired rank changes."""

import argparse
import itertools
import runpy
import sys
import time
from pathlib import Path

import pandas as pd
from scipy.stats import kendalltau, spearmanr

from run_transfer import ASSAY, BENCH, GROUPS

sys.path.insert(0, str(BENCH / "ligand_hydroxyl"))
from model import ARMS, Receptor, score

count = runpy.run_path(str(BENCH / "terminal_rotor_score/independent/verify.py"))[
    "independently_count"
]


def evaluate(work, out):
    out.mkdir(exist_ok=True, parents=True)
    geom = pd.read_csv(work / "geometry.csv")
    query = {
        s: (work / s / "ligand.pdbqt").read_text() for es in GROUPS.values() for s in es
    }
    scale = {s: 1 / (1 + 0.05846 * count(text.encode())) for s, text in query.items()}
    rec = {}
    times = []
    for r in query:
        start = time.perf_counter()
        rec[r] = Receptor((work / r / "receptor.pdbqt").read_text())
        times.append(
            dict(
                operation="build_receptor",
                receiver=r,
                query="",
                cap="",
                seed="",
                poses="",
                wall_s=time.perf_counter() - start,
            )
        )
    poses, cells = [], []
    for (q, r, cap, seed, engine), group in geom.groupby(
        ["query", "receiver", "cap", "seed", "engine"]
    ):
        key = dict(
            target=group.iloc[0].target,
            query=q,
            receiver=r,
            kind="self" if q == r else "cross",
            cap=int(cap),
            seed=int(seed),
            engine=engine,
        )
        path = work / "runs" / f"{q}_{r}_{cap}_{seed}_{engine}.pdbqt"
        blocks = [b for b in path.read_text().split("ENDMDL") if b.strip()]
        records = []
        start = time.perf_counter()
        for b, row in zip(blocks, group.sort_values("model").itertuples(), strict=True):
            rank = int(b.split("MODEL")[1].split()[0])
            if rank != row.model:
                raise ValueError("Wrong pose order")
            item = dict(
                **key,
                model=rank,
                rmsd=row.rmsd,
                native=row.dry,
                scale=scale[q],
                n_ca=row.n_ca,
                ca_rmsd=row.ca_rmsd,
            )
            if engine == "pro":
                terms = score(rec[r], query[q], b)
                item.update(
                    {a: row.dry + scale[q] * terms[a] for a in ARMS},
                    sites=terms["sites"],
                    active=terms["active_sites"],
                    contacts=terms["contacts"],
                )
            records.append(item)
        if engine == "pro":
            times.append(
                dict(
                    operation="score_pool",
                    receiver=r,
                    query=q,
                    cap=cap,
                    seed=seed,
                    poses=len(records),
                    wall_s=time.perf_counter() - start,
                )
            )
        poses.extend(records)
        for arm in ("native", *ARMS) if engine == "pro" else ("native",):
            # Native output order is authoritative when reports tie at 0.001.
            ranked = (
                sorted(records, key=lambda p: (p[arm], p["model"]))
                if arm != "native"
                else records
            )
            rms = [p["rmsd"] for p in ranked]
            cells.append(
                dict(
                    **key,
                    arm=arm,
                    score=ranked[0][arm],
                    index=ranked[0]["model"],
                    top1_A=rms[0],
                    top1=int(rms[0] <= 2),
                    top5=int(min(rms[:5]) <= 2),
                    oracle=int(min(rms) <= 2),
                    pool=len(rms),
                )
            )
    pd.DataFrame(poses).to_csv(out / "poses.csv", index=False)
    cells = pd.DataFrame(cells)
    cells.to_csv(out / "cells.csv", index=False)
    pd.DataFrame(times).to_csv(out / "rescore_cost.csv", index=False)
    summary = cells.groupby(
        ["target", "kind", "cap", "engine", "arm"], as_index=False
    ).agg(
        pools=("pool", "count"),
        top1=("top1", "sum"),
        top5=("top5", "sum"),
        oracle=("oracle", "sum"),
    )
    summary["conditional_top1"] = summary.top1 / summary.oracle
    summary.to_csv(out / "summary.csv", index=False)
    pair_scores = []
    affinities = []
    for (cap, engine, arm), frame in cells[
        (cells.target == "ERRg") & (cells.receiver == "6I63")
    ].groupby(["cap", "engine", "arm"]):
        pred = frame.groupby("query").score.median()
        for population, ids in [
            ("single_microstate", [s for s in ASSAY if s != "6I66"]),
            ("including_6I66", list(ASSAY)),
        ]:
            eligible = []
            for s, t in itertools.combinations(sorted(ids, key=ASSAY.get), 2):
                factor = ASSAY[t] / ASSAY[s]
                cred = 0.5 if pred[s] == pred[t] else float(pred[s] < pred[t])
                pair_scores.append(
                    dict(
                        population=population,
                        cap=cap,
                        engine=engine,
                        arm=arm,
                        stronger=s,
                        weaker=t,
                        measured_ratio=factor,
                        pred_stronger=pred[s],
                        pred_weaker=pred[t],
                        credit=cred,
                        eligible=int(factor >= 3),
                    )
                )
                if factor >= 3:
                    eligible.append(cred)
            energies = pred.loc[ids]
            observed = [ASSAY[s] for s in ids]
            affinities.append(
                dict(
                    population=population,
                    cap=cap,
                    engine=engine,
                    arm=arm,
                    n=len(ids),
                    eligible=len(eligible),
                    correct=sum(eligible),
                    spearman=spearmanr(observed, energies).statistic,
                    kendall=kendalltau(observed, energies).statistic,
                )
            )
    pd.DataFrame(pair_scores).to_csv(out / "affinity_pairs.csv", index=False)
    pd.DataFrame(affinities).to_csv(out / "affinity.csv", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    evaluate(a.work, a.output)
