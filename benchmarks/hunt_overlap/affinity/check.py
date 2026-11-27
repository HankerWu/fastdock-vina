"""Independent readback of all native/parent/candidate/control two-state consensus scores."""

import csv
import itertools
import json
import re
import sys
from collections import defaultdict
from decimal import Decimal
from pathlib import Path

repo, work = map(Path, sys.argv[1:])
with (repo / "benchmarks/galectin_affinity/table1.csv").open() as f:
    paper = list(csv.DictReader(f))
states = {"Gal1": ["1GZW", "4Q26"], "Gal3": ["1KJL", "3T1L"]}
with (work / "rankings.csv").open() as f:
    summary = list(csv.DictReader(f))
cases = defaultdict(dict)
for path in (work / "raw").glob("*.pdbqt"):
    rec, cid, seed, arm, cap = path.stem.split("-")
    score = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", path.read_text()).group(1)
    cases[f"{arm}-{seed}-{cap}"][rec, cid] = Decimal(score)
consensuses = 0
for row in summary:
    if (
        row["state"] != "consensus"
        or row["subset"] != "all"
        or row["panel"] == "heavy_atoms"
    ):
        continue
    target = row["target"]
    preds = cases[row["panel"]]
    total = 0
    correct = Decimal(0)
    for a, b in itertools.permutations(paper, 2):
        if Decimal(a[target + "_Kd_uM"]) * 3 > Decimal(b[target + "_Kd_uM"]):
            continue
        gap = sum(
            preds[rec, b["compound"]] - preds[rec, a["compound"]]
            for rec in states[target]
        )
        total += 1
        correct += Decimal(1) if gap > 0 else Decimal(".5") if gap == 0 else Decimal(0)
    assert (total, correct) == (int(row["pairs"]), Decimal(row["correct"])), row
    consensuses += 1
print(
    json.dumps(
        {
            "validated_consensus_counts": consensuses,
            "pdbqts": sum(map(len, cases.values())),
        }
    )
)
