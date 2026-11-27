"""Recompute every developmental pose and counter and tabulate the paired endpoints."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import sys


BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH / "rigid_flip"))
from holdout import paired

sys.path.insert(0, str(BENCH / "bond_distance/pose"))
from audit import independent_rmsds


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def analyze(work):
    totals = []
    details = []
    for study in ("primary", "reduced"):
        d = work / study
        rows = read(d / "cells.csv")
        expected = {
            (m["key"], int(m["mode"])): float(m["rmsd_A"])
            for m in read(d / "modes.csv")
        }
        for r in rows:
            q, rec = r["query"], r["receiver"]
            meta = {
                s: json.loads((d / s / "preparation.json").read_text())
                for s in (q, rec)
            }
            stem = d / "runs" / r["key"]
            output = stem.with_suffix(".pdbqt")
            assert hashlib.sha256(output.read_bytes()).hexdigest() == r["sha256"]
            values = independent_rmsds(meta[q], meta[rec], output)
            assert len(values) == int(r["modes"])
            assert (
                max(abs(v - expected[r["key"], i]) for i, v in enumerate(values, 1))
                < 1e-7
            )
            log = stem.with_suffix(
                ".counts" if r["arm"] == "upstream" else ".log"
            ).read_text()
            assert all(int(r[k]) == v for k, v in paired.counts(log).items())
            measures = dict(re.findall(r"\b(axes|other|offered|redirected)=(\d+)", log))
            details.append(
                dict(
                    study=study,
                    key=r["key"],
                    **{
                        k: measures.get(k, "")
                        for k in ("axes", "other", "offered", "redirected")
                    },
                )
            )
        for target in ("all", *sorted({r["target"] for r in rows})):
            for arm in sorted({r["arm"] for r in rows}):
                group = [
                    r
                    for r in rows
                    if r["arm"] == arm and target in ("all", r["target"])
                ]
                totals.append(
                    dict(
                        study=study,
                        target=target,
                        arm=arm,
                        n=len(group),
                        **{
                            k: sum(float(r[k]) <= 2 for r in group)
                            for k in ("top1_A", "best5_A", "oracle_A")
                        },
                        clashes=sum(int(r["top1_clash_pairs"]) > 0 for r in group),
                        evals=sum(
                            int(r[k])
                            for r in group
                            for k in ("hunt_evals", "second_evals", "explicit_evals")
                        ),
                        **{
                            k: sum(float(r[k]) for r in group)
                            for k in ("wall_s", "cpu_s")
                        },
                        max_rss_KiB=max(int(r["peak_rss_KiB"]) for r in group),
                    )
                )
    paired.write(work / "structure_summary.csv", totals)
    paired.write(work / "torsion_counts.csv", details)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    analyze(parser.parse_args().work)
