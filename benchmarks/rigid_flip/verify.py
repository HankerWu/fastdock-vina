"""Independently reproduce every archived pose RMSD and check primary search counts and hits."""

import csv
import hashlib
import json
import re
import sys
import tarfile
import tempfile
from pathlib import Path

from holdout import BENCH, CASES

sys.path.insert(0, str(BENCH / "bond_distance/pose"))
from audit import independent_rmsds


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    report = {}
    with tempfile.TemporaryDirectory() as tmp:
        w = Path(tmp)
        with tarfile.open(here / "inputs.tar.xz") as archive:
            archive.extractall(w, filter="data")
        meta = {q: json.loads((w / q / "preparation.json").read_text()) for q in CASES}
        for site in json.loads((here / "inputs.json").read_text()):
            assert all(
                hashlib.sha256((w / name).read_bytes()).hexdigest() == sha
                for name, sha in site["files"].items()
            )
        for study, expected in (("primary", 112), ("reduced", 56), ("newseeds", 56)):
            source = here / study
            dest = w / study
            with tarfile.open(source / "raw.tar.xz") as archive:
                archive.extractall(dest, filter="data")
            rows = read(source / "cells.csv")
            modes = {
                (m["key"], int(m["mode"])): float(m["rmsd_A"])
                for m in read(source / "modes.csv")
            }
            assert len(rows) == expected and len({r["key"] for r in rows}) == expected
            observed = 0
            counts = {}
            worst = 0
            for r in rows:
                k = r["key"]
                prefix = dest / "runs" / k
                pdbqt = prefix.with_suffix(".pdbqt")
                assert r["sha256"] == hashlib.sha256(pdbqt.read_bytes()).hexdigest()
                vals = independent_rmsds(meta[r["query"]], meta[r["receiver"]], pdbqt)
                observed += len(vals)
                worst = max(
                    [worst] + [abs(v - modes[k, i]) for i, v in enumerate(vals, 1)]
                )
                assert all(
                    abs(float(r[field]) - v) < 1e-7
                    for field, v in zip(
                        ("top1_A", "best5_A", "oracle_A"),
                        (vals[0], min(vals[:5]), min(vals)),
                    )
                )
                log = prefix.with_suffix(
                    ".counts" if r["arm"] == "upstream" else ".log"
                ).read_text()
                items = dict(
                    re.findall(
                        r"\b((?:hunt|second|explicit)_(?:calls|evals))=(\d+)", log
                    )
                )
                assert len(items) == 6 and all(
                    r[name] == value for name, value in items.items()
                ), k
                hits = counts.setdefault(r["arm"], [0, 0, 0])
                for i, v in enumerate((vals[0], min(vals[:5]), min(vals))):
                    hits[i] += v <= 2
            assert len(modes) == observed and worst < 1e-7
            report[study] = {
                "runs": len(rows),
                "poses": observed,
                "max_rmsd_delta_A": worst,
                "hits": counts,
            }
    print(json.dumps(report, indent=2))
