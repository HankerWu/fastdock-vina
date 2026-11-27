"""Pair each lower-cap result with the same parent's target and random seed."""

import argparse
from pathlib import Path

from external_audit import rows
from external_run import write


def cost(low, base, base_repeat, out, low_repeat=None):
    repeated = {r["key"]: r for r in rows(base_repeat / "runs.csv")}
    parent = {
        (r["pdb"], r["seed"], r["cap"]): r
        for r in rows(base / "structures.csv")
        if r["arm"] == "parent"
    }
    repeats = {r["key"]: r for r in rows(low_repeat / "runs.csv")} if low_repeat else {}
    pairs = []
    for r in rows(low / "structures.csv"):
        if r["cap"] not in ("9000", "2000"):
            continue
        ref = parent[r["pdb"], r["seed"], {"2000": "3000", "9000": "12000"}[r["cap"]]]
        base_ns = (int(ref["elapsed_ns"]) + int(repeated[ref["key"]]["elapsed_ns"])) / 2
        low_ns = (
            (int(r["elapsed_ns"]) + int(repeats[r["key"]]["elapsed_ns"])) / 2
            if repeats
            else int(r["elapsed_ns"])
        )
        pairs.append(
            {
                "key": r["key"],
                "parent_key": ref["key"],
                "eval_ratio": round(int(r["joint"]) / int(ref["joint"]), 6),
                "wall_ratio": round(low_ns / base_ns, 6),
                "top1_A": r["top1_A"],
                "parent_top1_A": ref["top1_A"],
                "oracle_A": r["oracle_A"],
                "parent_oracle_A": ref["oracle_A"],
            }
        )
    write(out, pairs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("low", "parent", "parent_repeats", "out"):
        p.add_argument(k, type=Path)
    p.add_argument("--low-repeat", type=Path)
    a = p.parse_args()
    cost(a.low, a.parent, a.parent_repeats, a.out, a.low_repeat)
