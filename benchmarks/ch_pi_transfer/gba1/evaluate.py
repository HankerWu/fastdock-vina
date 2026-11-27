"""Rescore the completed searches and evaluate symmetry-corrected poses in receptor coordinates."""

import argparse
import csv
import itertools
import json
import runpy
import sys
import time
from pathlib import Path

import numpy as np
from experiment import CASES, HERE, cross, write_csv

sys.path.insert(0, str(HERE.parent))
from model import CHPi, WEIGHT

count_torsions = runpy.run_path(
    str(HERE.parents[1] / "terminal_rotor_score/independent/verify.py")
)["independently_count"]


def distribution(values):
    q = np.quantile(values, [0.25, 0.5, 0.75, 0.90, 0.95])
    return dict(
        mean=float(np.mean(values)),
        q25=q[0],
        median=q[1],
        q75=q[2],
        p90=q[3],
        p95=q[4],
        max=max(values),
    )


def evaluate(w, out, cap=None):
    run_dir = w / f"budget_{cap}" if cap else w
    expected = 27 if cap else 81
    meta = json.loads((w / "metadata.json").read_text())
    with (run_dir / "compute.csv").open() as f:
        runs = list(csv.DictReader(f))
    if len(runs) != expected or len({r["args"] for r in runs}) != expected:
        raise ValueError("Expected a complete ligand/receptor/seed manifest")
    models = {}
    modes, cells, costs, geom = [], [], [], []
    for r in runs:
        src, dst, state, engine, seed = (
            r[k] for k in ("src", "dst", "state", "arm", "seed")
        )
        fields = dict(
            query=src,
            receptor=dst,
            kind="redock" if src == dst else "cross",
            state=state,
            seed=seed,
        )
        key = (src, dst, state)
        if key not in models:
            rx = (w / dst / "receptor.pdbqt").read_text()
            lx = (w / src / (state + ".pdbqt")).read_text()
            start = time.perf_counter()
            models[key] = CHPi(rx, lx), 1 / (1 + 0.05846 * count_torsions(lx.encode()))
            costs.append(
                dict(
                    **fields,
                    stage="instantiate",
                    n=1,
                    wall_s=time.perf_counter() - start,
                )
            )
        model, factor = models[key]
        transform, shift, residual, n = cross.alignment(meta[src], meta[dst])
        if not any(g["query"] == src and g["receptor"] == dst for g in geom):
            geom.append(
                dict(
                    query=src, receptor=dst, pocket_ca_count=n, pocket_fit_rmsd=residual
                )
            )
        qmeta = meta[src]
        refs = np.array(qmeta["mapped"]) @ transform + shift
        labels = qmeta["serial" if state == "neutral" else "protonated_serial"]
        ranks = ("native", "chpi") if engine == "pro" else ("official",)
        path = run_dir / "runs" / (f"{src}_{dst}_{state}_{engine}_{seed}.pdbqt")
        pool = []
        blocks = path.read_text().split("ENDMDL") if r["status"] == "ok" else []
        elapsed = 0.0
        for rank, block in enumerate((b for b in blocks if b.strip()), 1):
            xyz = {
                int(l[6:11]): cross.coords(l)
                for l in block.splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
            }
            points = np.array([xyz[k] for k in labels])
            errors = np.sqrt(np.mean(np.sum((points - refs) ** 2, axis=2), axis=1))
            native = float(
                next(
                    l.split()[3]
                    for l in block.splitlines()
                    if l.startswith("REMARK VINA RESULT:")
                )
            )
            corrected = forward = reverse = ""
            if engine == "pro":
                start = time.perf_counter()
                forward, reverse = model.components(block)
                elapsed += time.perf_counter() - start
                corrected = native + WEIGHT * factor * (forward + reverse)
            row = dict(
                **fields,
                engine=engine,
                rank=rank,
                rmsd=float(min(errors)),
                native=native,
                chpi=corrected,
                raw_forward=forward,
                raw_reverse=reverse,
                rotor_scale=factor,
            )
            pool.append(row)
            modes.append(row)
        if engine == "pro":
            costs.append(
                dict(**fields, stage="score_pool", n=len(pool), wall_s=elapsed)
            )
        for arm in ranks:
            entry = dict(
                **fields,
                arm=arm,
                count=len(pool),
                original_rank="",
                top1="",
                best5="",
                pool_min="",
                regret="",
            )
            if pool:
                ordered = (
                    sorted(pool, key=lambda p: p["chpi"]) if arm == "chpi" else pool
                )
                top = ordered[0]
                oracle = min(p["rmsd"] for p in pool)
                entry.update(
                    original_rank=top["rank"],
                    top1=top["rmsd"],
                    best5=min(p["rmsd"] for p in ordered[:5]),
                    pool_min=oracle,
                    regret=top["rmsd"] - oracle,
                )
            cells.append(entry)
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (
        ("modes", modes),
        ("cells", cells),
        ("rerank_cost", costs),
        ("alignment", geom),
    ):
        write_csv(out / (name + ".csv"), data)

    sums, pairs, ecdf, boot = [], [], [], []
    for state, subset in itertools.product(
        dict.fromkeys(c["state"] for c in cells), ("all", "redock", "cross", *CASES)
    ):
        selected = [
            c
            for c in cells
            if c["state"] == state and subset in ("all", c["kind"], c["query"])
        ]
        byarm = {
            arm: [c for c in selected if c["arm"] == arm]
            for arm in sorted({c["arm"] for c in selected})
        }
        for arm, rows in byarm.items():
            valid = [c for c in rows if c["count"]]
            if not valid:
                continue
            for metric in ("top1", "best5", "pool_min", "regret"):
                values = [c[metric] for c in valid]
                sums.append(
                    dict(
                        state=state,
                        subset=subset,
                        arm=arm,
                        metric=metric,
                        n=len(values),
                        failures=len(rows) - len(values),
                        **distribution(values),
                    )
                )
            values = sorted(c["top1"] for c in valid)
            ecdf.extend(
                dict(
                    state=state,
                    subset=subset,
                    arm=arm,
                    order=i,
                    rmsd=v,
                    fraction=i / len(values),
                )
                for i, v in enumerate(values, 1)
            )
        changes = []
        for old, new in zip(byarm["native"], byarm["chpi"], strict=True):
            assert all(old[k] == new[k] for k in ("query", "receptor", "seed", "count"))
            if old["count"]:
                delta = new["top1"] - old["top1"]
                changes.append((old["seed"], delta))
                if subset == "all":
                    pairs.append(
                        dict(
                            state=state,
                            query=old["query"],
                            receptor=old["receptor"],
                            seed=old["seed"],
                            native=old["top1"],
                            chpi=new["top1"],
                            delta=delta,
                            chosen_rank=new["original_rank"],
                        )
                    )
        if not changes:
            continue
        deltas = [d for _, d in changes]
        seeds = sorted({s for s, _ in changes})
        effects = [
            np.mean([d for s in draw for obs, d in changes if s == obs])
            for draw in itertools.product(seeds, repeat=len(seeds))
        ]
        ci = np.quantile(effects, [0.025, 0.975])
        boot.append(
            dict(
                state=state,
                subset=subset,
                improved=sum(d < 0 for d in deltas),
                worsened=sum(d > 0 for d in deltas),
                same=sum(d == 0 for d in deltas),
                ci025=ci[0],
                ci975=ci[1],
                **distribution(deltas),
            )
        )
    for name, data in (
        ("distributions", sums),
        ("paired", pairs),
        ("ecdf", ecdf),
        ("paired_summary", boot),
    ):
        write_csv(out / (name + ".csv"), data)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--cap", type=int, choices=(12100, 24000))
    args = p.parse_args()
    evaluate(args.work.resolve(), args.output.resolve(), args.cap)
