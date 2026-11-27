"""Interleave original and reduced caps on one pinned CPU; assert primary coordinates."""

import argparse
import csv
import hashlib
import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

from holdout import CASES, args_for, paired

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("primary", type=Path)
    ap.add_argument("reduced", type=Path)
    for n in ("parent", "always", "adaptive"):
        ap.add_argument("--" + n, type=Path, required=True)
    a = ap.parse_args()
    d = a.primary.resolve()
    meta = {q: json.loads((d / q / "preparation.json").read_text()) for q in CASES}
    with (d / "cells.csv").open() as f:
        cells = [r for r in csv.DictReader(f) if r["arm"] == "parent"]
    schedule = random.Random(324901)
    result = []
    with tempfile.TemporaryDirectory() as raw:
        scratch = Path(raw)
        log, stats, out = [scratch / p for p in ("log", "cost", "output")]
        for rep in range(3):
            schedule.shuffle(cells)
            for cell in cells:
                choices = [
                    ("primary", "parent"),
                    ("primary", "always"),
                    ("primary", "adaptive"),
                    ("reduced", "always"),
                    ("reduced", "adaptive"),
                ]
                schedule.shuffle(choices)
                for design, arm in choices:
                    cap, seed = int(cell["cap"]), int(cell["seed"])
                    w = d if design == "primary" else a.reduced
                    q, r = cell["query"], cell["receiver"]
                    budget = cap if design == "primary" else cap * 85 // 100
                    cli = args_for(d, meta, q, r, seed, budget, out)
                    subprocess.run(
                        [
                            sys.executable,
                            "-S",
                            str(Path(__file__).with_name("process_cost.py")),
                            str(log),
                            str(stats),
                            str(getattr(a, arm).resolve()),
                            *cli,
                            "--search_counters",
                        ],
                        check=True,
                    )
                    key = f"{q}_{r}_{seed}_{budget}_{arm}"
                    original = json.loads((w / "runs" / (key + ".json")).read_text())
                    if (
                        hashlib.sha256(out.read_bytes()).hexdigest()
                        != original["sha256"]
                    ):
                        raise ValueError(key)
                    result.append(
                        {
                            "rep": rep,
                            "query": q,
                            "receiver": r,
                            "seed": seed,
                            "parent_cap": cap,
                            "cap": budget,
                            "arm": arm,
                            "design": design,
                            **json.loads(stats.read_text()),
                            "sha256": original["sha256"],
                        }
                    )
            paired.write(d / "cost_replay.csv", result)
