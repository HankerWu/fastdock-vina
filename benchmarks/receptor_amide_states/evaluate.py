"""Score all development pools in the preregistered archived transfer panel."""

import argparse
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
import pandas as pd
from model import RT, AmideReceptor

HERE = Path(__file__).resolve().parent
KEY = ["target", "query", "receiver", "kind", "cap", "seed"]
ARMS = ["native", "flip_all", "minimum", "uniform", "ensemble"]


def quantiles(values):
    out = {"mean": values.mean()}
    out.update(
        zip(
            ("q25", "median", "q75", "p90", "p95"),
            np.quantile(values, [0.25, 0.5, 0.75, 0.90, 0.95]),
        )
    )
    return out


def evaluate(work, out):
    dataset = pd.read_csv(HERE.parent / "ligand_hydroxyl/transfer/results/poses.csv")
    dataset = dataset[dataset.engine == "pro"]
    out.mkdir(exist_ok=True, parents=True)
    recs, poses, cells, sites, weights, costs = {}, [], [], [], [], []
    for key, group in dataset.groupby(KEY):
        info = dict(zip(KEY, key))
        _, q, r, _, cap, seed = key
        if r not in recs:
            wall, cpu = perf_counter(), process_time()
            rec = AmideReceptor((work / r / "receptor.pdbqt").read_text())
            recs[r] = rec
            costs.append(
                {
                    "receiver": r,
                    "action": "receptor",
                    "sites": len(rec.sites),
                    "excluded": len(rec.excluded),
                    "wall_s": perf_counter() - wall,
                    "cpu_s": process_time() - cpu,
                }
            )
            for s in rec.sites:
                weights.append(
                    {
                        "receiver": r,
                        "residue": s.residue,
                        "apo_delta": s.apo_delta,
                        "p_flip": 1 / (1 + np.exp(s.apo_delta / RT)),
                    }
                )
        text = (work / q / "ligand.pdbqt").read_text()
        blocks = [
            b
            for b in (work / "runs" / f"{q}_{r}_{cap}_{seed}_pro.pdbqt")
            .read_text()
            .split("ENDMDL")
            if b.strip()
        ]
        wall, cpu = perf_counter(), process_time()
        records = []
        for row, block in zip(
            group.sort_values("model").itertuples(), blocks, strict=True
        ):
            native = float(
                next(
                    l for l in block.splitlines() if l.startswith("REMARK VINA RESULT:")
                ).split()[3]
            )
            if native != row.native or row.model != int(
                block.split("MODEL")[1].split()[0]
            ):
                raise ValueError("Archived score/pose identity mismatch")
            scores = recs[r].score(text, block, detail=True)
            for resid, apo, dl, correction in scores.pop("sites"):
                sites.append(
                    dict(
                        **info,
                        model=row.model,
                        residue=resid,
                        apo_delta=apo,
                        ligand_delta=dl,
                        correction=correction,
                    )
                )
            record = dict(
                **info,
                model=row.model,
                rmsd_A=row.rmsd,
                scale=row.scale,
                native=native,
                active_sites=scores["active_sites"],
            )
            for arm in ARMS[1:]:
                record[arm] = native + row.scale * scores[arm]
            records.append(record)
        costs.append(
            dict(
                **info,
                action="pool",
                poses=len(records),
                wall_s=perf_counter() - wall,
                cpu_s=process_time() - cpu,
            )
        )
        poses.extend(records)
        for arm in ARMS:
            sorted_rows = (
                records
                if arm == "native"
                else sorted(records, key=lambda rr: (rr[arm], rr["model"]))
            )
            rms = [p["rmsd_A"] for p in sorted_rows]
            cells.append(
                dict(
                    **info,
                    arm=arm,
                    count=len(rms),
                    selected=sorted_rows[0]["model"],
                    score=sorted_rows[0][arm],
                    top1_A=rms[0],
                    best5_A=min(rms[:5]),
                    oracle_A=min(rms),
                    regret_A=rms[0] - min(rms),
                    delta_A=rms[0] - records[0]["rmsd_A"],
                )
            )
    for name, data in [
        ("poses", poses),
        ("cells", cells),
        ("site_changes", sites),
        ("priors", weights),
        ("cost", costs),
    ]:
        pd.DataFrame(data).to_csv(out / f"{name}.csv", index=False)
    data = pd.DataFrame(cells)
    summary = []
    for target, sub in [("all", data), *data.groupby("target")]:
        for arm, frame in sub.groupby("arm"):
            d = frame.delta_A
            for field in ["top1_A", "best5_A", "oracle_A", "regret_A", "delta_A"]:
                summary.append(
                    dict(
                        target=target,
                        arm=arm,
                        field=field,
                        n=len(frame),
                        better=(d < 0).sum(),
                        worse=(d > 0).sum(),
                        unchanged=(d == 0).sum(),
                        **quantiles(frame[field]),
                    )
                )
    pd.DataFrame(summary).to_csv(out / "summary.csv", index=False)
    # Resample complete ligand clusters within family; preserve both seeds, budgets and receptors.
    rng = np.random.default_rng(191009)
    cis = []
    for target, sub in data[data.arm == "ensemble"].groupby("target"):
        clusters = [g.delta_A.to_numpy() for _, g in sub.groupby("query")]
        samples = [
            np.mean(np.concatenate([clusters[j] for j in draw]))
            for draw in rng.integers(0, len(clusters), (5000, len(clusters)))
        ]
        lo, hi = np.quantile(samples, [0.025, 0.975])
        cis.append(
            {
                "target": target,
                "queries": len(clusters),
                "paired_mean_A": sub.delta_A.mean(),
                "q025_A": lo,
                "q975_A": hi,
            }
        )
    pd.DataFrame(cis).to_csv(out / "query_cluster_interval.csv", index=False)
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    evaluate(args.work, args.out)
