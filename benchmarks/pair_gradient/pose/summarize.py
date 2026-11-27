"""Reduce paired structural cells without treating seed repeats as new complexes."""

import argparse
from pathlib import Path

import pandas as pd

COUNTS = ["hunt_evals", "second_evals", "explicit_evals"]


def summarize(work):
    rows = pd.concat(
        [
            pd.read_csv(work / f)
            for f in ("cells.csv", "shift_cells.csv", "release_cells.csv")
        ]
    )
    rows["joint"] = rows[COUNTS].sum(axis=1, min_count=3)
    rows["local_calls"] = rows[["hunt_calls", "second_calls", "explicit_calls"]].sum(
        axis=1, min_count=3
    )
    assert rows[rows.arm != "release"].joint.notna().all()
    parent = rows[rows.arm == "parent"].set_index(["query", "receiver", "seed", "cap"])
    summaries = []
    for keys in (
        ["arm"],
        ["cap", "arm"],
        ["target", "arm"],
        ["target", "kind", "arm"],
        ["query", "receiver", "arm"],
    ):
        for _, g in rows.groupby(keys):
            d = {name: g[name].iloc[0] for name in keys}
            p = parent.loc[
                list(zip(g["query"], g.receiver, g.seed, g.cap))
            ].reset_index(drop=True)
            g = g.reset_index(drop=True)
            hit, baseline = g.top1_A <= 2, p.top1_A <= 2
            qualified = (g.seconds <= p.seconds) & (g.joint <= p.joint)
            summaries.append(
                dict(
                    **d,
                    n=len(g),
                    top1=hit.sum(),
                    top5=(g.best5_A <= 2).sum(),
                    oracle=(g.oracle_A <= 2).sum(),
                    gained=(hit & ~baseline).sum(),
                    lost=(~hit & baseline).sum(),
                    top1_A_median=g.top1_A.median(),
                    oracle_A_median=g.oracle_A.median(),
                    joint=g.joint.sum(min_count=1),
                    local_calls=g.local_calls.sum(min_count=1),
                    wall_s=g.seconds.sum(),
                    top1_with_clash=(g.top1_clash_pairs > 0).sum(),
                    modes=g.modes.sum(),
                    lower_time_and_joint=qualified.sum(),
                    gains_in_lower_cost=(hit & ~baseline & qualified).sum(),
                )
            )
    pd.DataFrame(summaries).to_csv(work / "summary.csv", index=False)
    matched = pd.read_csv(work / "matched_cells.csv")
    matched["joint"] = matched[COUNTS].sum(axis=1)
    control = parent.reset_index()
    control = control[(control.seed == 98251) & (control.cap == 12000)]
    merged = matched.merge(
        control[["query", "receiver", "key", "seconds", "joint", "top1_A"]],
        on=["query", "receiver"],
        suffixes=("", "_parent"),
        validate="many_to_one",
    )
    repeat = (
        pd.read_csv(work / "time_replay.csv")
        .groupby("key")[["seconds", "cpu_seconds"]]
        .median()
    )
    for k in ("seconds", "cpu_seconds"):
        merged["replay_" + k] = merged.key.map(repeat[k])
        merged["replay_parent_" + k] = merged.key_parent.map(repeat[k])
    merged.to_csv(work / "matched_pairs.csv", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    summarize(p.parse_args().work)
