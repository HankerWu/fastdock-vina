"""Rebuild matched pair credits directly from output remarks and frozen primary Kd labels."""

import argparse
import csv
import itertools
import json
from pathlib import Path

KD_NM = {"4OGI": 37, "4O76": 123, "4OGJ": 164, "4O77": 3400}


def check(work, controls, results):
    with (results / "pairs.csv").open() as f:
        rows = {
            (r["receptor"], r["seed"], r["arm"], r["strong"], r["weak"]): r
            for r in csv.DictReader(f)
        }
    arms = ("native", "multi", "contact")
    hits = {a: 0.0 for a in arms}
    counts = {a: 0 for a in arms}
    for r, seed, arm in itertools.product(
        ("3MXF", "4OGI", "4OGJ"), ("314159", "271828", "161803", "141421"), arms
    ):
        score = {}
        for ligand in KD_NM:
            output = (
                work if arm == "contact" else controls
            ) / f"{r}_{ligand}_{seed}_{arm}.pdbqt"
            score[ligand] = next(
                float(l.split(":", 1)[1].split()[0])
                for l in output.read_text().splitlines()
                if l.startswith("REMARK VINA RESULT")
            )
        for strong, weak in itertools.permutations(KD_NM, 2):
            if KD_NM[weak] >= 3 * KD_NM[strong]:
                row = rows.pop((r, seed, arm, strong, weak))
                margin = score[weak] - score[strong]
                credit = (margin > 0) + 0.5 * (margin == 0)
                if abs(margin - float(row["score_margin"])) > 1e-9 or credit != float(
                    row["credit"]
                ):
                    raise ValueError("Incorrect margin/credit")
                hits[arm] += credit
                counts[arm] += 1
    if rows:
        raise ValueError("Unexpected pair rows")
    report = {"pair_counts": counts, "pair_credits": hits}
    print(json.dumps(report))
    (results / "pair_check.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("controls", type=Path)
    p.add_argument("results", type=Path)
    a = p.parse_args()
    check(a.work, a.controls, a.results)
