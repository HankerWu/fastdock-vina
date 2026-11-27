"""Replay the first switched proposal in selected post-hoc trajectories, with identical output checks."""

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

from holdout import BENCH, CASES, args_for, paired

sys.path.insert(0, str(BENCH / "bond_distance/pose"))
from audit import independent_rmsds

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument("--always", type=Path, required=True)
    ap.add_argument("--adaptive", type=Path, required=True)
    a = ap.parse_args()
    w = a.work.resolve()
    meta = {q: json.loads((w / q / "preparation.json").read_text()) for q in CASES}
    out = w / "first_branch"
    out.mkdir(exist_ok=True)
    records = []
    for q, r, s, c in (
        ("3ERT", "1SJ0", 202711, 3000),
        ("3ERT", "1SJ0", 202711, 12000),
        ("4UND", "4UND", 202727, 12000),
        ("3ERT", "1SJ0", 202727, 12000),
        ("1T49", "1T48", 202711, 12000),
    ):
        for arm in ("adaptive", "always"):
            key = f"{q}_{r}_{s}_{c}_{arm}"
            prefix = out / key
            probe = out / (key + ".output")
            log = out / (key + ".log")
            with log.open("w") as f:
                subprocess.run(
                    [
                        str(getattr(a, arm).resolve()),
                        *args_for(w, meta, q, r, s, c, probe),
                        "--search_counters",
                    ],
                    check=True,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    env=dict(os.environ, SHAPE_TRACE=str(prefix)),
                )
            expected = (w / "runs" / (key + ".pdbqt")).read_bytes()
            if probe.read_bytes() != expected:
                raise ValueError("Observer changed output: " + key)
            digest = hashlib.sha256(expected).hexdigest()
            probe.unlink()
            raw = prefix.with_suffix(".pdbqt")
            if not raw.exists():
                records.append(
                    {
                        "key": key,
                        "branch": "unswitched",
                        "step": "",
                        "source_A": "",
                        "offer_A": "",
                        "hunt_A": "",
                        "accepted": "",
                        "evals": "",
                        "previous_energy": "",
                        "hunt_energy": "",
                        "sha256": digest,
                    }
                )
                continue
            rms = independent_rmsds(meta[q], meta[r], raw)
            assert len(rms) == 5
            with prefix.with_suffix(".tsv").open() as f:
                for i, row in enumerate(csv.DictReader(f, delimiter="\t")):
                    records.append(
                        {
                            "key": key,
                            "branch": row["branch"],
                            "step": row["step"],
                            "source_A": rms[0],
                            "offer_A": rms[1 + i * 2],
                            "hunt_A": rms[2 + i * 2],
                            "accepted": row["accepted"],
                            "evals": row["evals"],
                            "previous_energy": row["previous_energy"],
                            "hunt_energy": row["hunt_energy"],
                            "sha256": digest,
                        }
                    )
    paired.write(w / "first_branch.csv", records)
