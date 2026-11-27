"""Check deterministic controls and report paired receptor and graph effects, in Å."""

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from summarize import read, save


def evaluate(run, dest):
    dest.mkdir(exist_ok=True, parents=True)
    pairs, distributions, controls = [], [], []
    historic = {
        (r["arm"], r["drug"], r["seed"]): r
        for r in read(HERE.parent / "results/cells.csv")
    }
    for cohort in ("historical", "fresh"):
        rows = read(run / cohort / "cells.csv")
        by_key = {(r["arm"], r["drug"], r["seed"]): r for r in rows}
        if len(rows) != len(by_key):
            raise ValueError("Duplicate docking cell")
        if cohort == "historical":
            for (arm, drug, seed), row in historic.items():
                if arm in ("official", "double"):
                    by_key[arm, drug, seed] = row
        drugs = sorted({x["drug"] for x in rows})
        seeds = sorted({x["seed"] for x in rows})
        for drug, seed in itertools.product(drugs, seeds):

            def check(name, x, y):
                controls.append(
                    {
                        "cohort": cohort,
                        "drug": drug,
                        "seed": seed,
                        "comparison": name,
                        "same_pool_bytes": int(x["pool_sha256"] == y["pool_sha256"]),
                        "same_evals": int(x["evals"] == y["evals"]),
                    }
                )

            def cell(arm):
                return by_key[arm, drug, seed]

            check("short_only vs far", cell("short_only"), cell("far_both"))
            check("long_only vs both", cell("long_only"), cell("both"))
            for arm in ("far_both", "both"):
                check(arm + " graph vs inferred", cell(arm + "_graph"), cell(arm))
            if cohort == "historical":
                for arm, ref in [("far_both", "dry"), ("both", "1OL7_released")]:
                    check(arm + " vs archived", cell(arm), historic[ref, drug, seed])
            for baseline, variant in (
                ("far_both", "both"),
                ("far_both_graph", "both_graph"),
                ("far_both", "far_both_graph"),
                ("both", "both_graph"),
            ):
                left, right = cell(baseline), cell(variant)
                pairs.append(
                    dict(
                        cohort=cohort,
                        drug=drug,
                        seed=seed,
                        baseline=baseline,
                        variant=variant,
                        **{
                            f"delta_{m}_A": float(right[m]) - float(left[m])
                            for m in ("top1", "best5", "oracle", "regret")
                        },
                    )
                )
        for arm in sorted({k[0] for k in by_key}):
            for drug in ("all", *drugs):
                group = [
                    r
                    for (a, d, _), r in by_key.items()
                    if a == arm and drug in ("all", d)
                ]
                errors = np.array([float(r["top1"]) for r in group])
                q1, median, q3, p90, p95 = np.percentile(errors, [25, 50, 75, 90, 95])
                distributions.append(
                    dict(
                        cohort=cohort,
                        arm=arm,
                        drug=drug,
                        n=len(group),
                        mean=errors.mean(),
                        q25=q1,
                        median=median,
                        q75=q3,
                        p90=p90,
                        p95=p95,
                        **{
                            m + "_mean": np.mean([float(r[m]) for r in group])
                            for m in ("best5", "oracle", "regret")
                        },
                        severe_ge_8=int((errors >= 8).sum()),
                        above_10=int((errors > 10).sum()),
                        empty=sum(int(r["n"]) == 0 for r in group),
                        modes=sum(int(r["n"]) for r in group),
                    )
                )
    save(dest / "controls.csv", controls)
    save(dest / "pairs.csv", pairs)
    save(dest / "distributions.csv", distributions)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    evaluate(a.run, a.out)
