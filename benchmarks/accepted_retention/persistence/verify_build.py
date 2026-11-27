"""Check a fresh recurrence build against every stored coordinate file and count."""

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "coverage"))
import tarfile

from probe import command, counters

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("archived", "structures", "results", "binary", "output"):
        p.add_argument(name, type=Path)
    p.add_argument("--cpu", type=int, required=True)
    a = p.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    count = 0
    records = csv.DictReader((a.results / "cells.csv").open())
    with (
        tempfile.TemporaryDirectory() as tmp,
        tarfile.open(a.results / "runs.tar.xz") as t,
    ):
        for r in records:
            if r["arm"] != "recurrence":
                continue
            w = a.structures if r["target"] in ("NUDT1", "CHIT1") else a.archived
            out = Path(tmp) / "pose.pdbqt"
            result = subprocess.run(
                command(w, r, a.binary, out), check=True, capture_output=True
            )
            pred = t.extractfile(r["key"] + ".pdbqt").read()
            log = t.extractfile(r["key"] + ".log").read()
            if out.read_bytes() != pred or counters(log) != counters(result.stdout):
                raise ValueError(r["key"])
            count += 1
    assert count == 96
    a.output.write_text(
        json.dumps(
            {
                "cases": count,
                "identical_pdbqts": count,
                "identical_counters": count,
                "rebuilt_binary_sha256": hashlib.sha256(
                    a.binary.read_bytes()
                ).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )
