"""Reparse all modes and counters; aggregate paired decisions without pooling targets."""

import argparse
import csv
import hashlib
import re
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "coverage"))
from probe import counters, mode_distances, references, table


def main(a):
    rows = list(csv.DictReader((a.results / "cells.csv").open()))
    if len(rows) != 288 or len({r["key"] for r in rows}) != 288:
        raise ValueError("Expected 96 cells and three arms")
    with tarfile.open(a.results / "runs.tar.xz") as t:
        for r in rows:
            work = a.structures if r["target"] in ("NUDT1", "CHIT1") else a.archived
            serial, refs = references(work, r)
            data = t.extractfile(r["key"] + ".pdbqt").read()
            log = t.extractfile(r["key"] + ".log").read()
            scores = mode_distances(data, serial, refs)
            counts = dict(re.findall(rb"(\w+)=(\d+)", b" ".join(counters(log))))
            chk = [scores[0], min(scores[:5]), min(scores)]
            if any(
                abs(x - float(r[k])) > 1e-10
                for x, k in zip(chk, ("top1", "top5", "oracle"))
            ):
                raise ValueError("Raw mode audit failed")
            if hashlib.sha256(data).hexdigest() != r["sha256"] or sum(
                int(counts[k])
                for k in (b"hunt_evals", b"second_evals", b"explicit_evals")
            ) != int(r["joint"]):
                raise ValueError("Raw hash or evaluations differ")
    grouped = []
    for split in ("development", "new_seed"):
        selected = [r for r in rows if r["split"] == split]
        for target in ("ALL", *sorted({r["target"] for r in selected})):
            for arm in ("parent", "energy", "recurrence"):
                block = [
                    r
                    for r in selected
                    if r["arm"] == arm and (target == "ALL" or target == r["target"])
                ]
                grouped.append(
                    {
                        "split": split,
                        "target": target,
                        "arm": arm,
                        "n": len(block),
                        "top1": sum(float(r["top1"]) <= 2 for r in block),
                        "top5": sum(float(r["top5"]) <= 2 for r in block),
                        "oracle": sum(float(r["oracle"]) <= 2 for r in block),
                        "joint": sum(int(r["joint"]) for r in block),
                        "local_calls": sum(int(r["local_calls"]) for r in block),
                        "slots": sum(int(r["slots"]) for r in block),
                        "replacements": sum(int(r["replaced"]) for r in block),
                        "wall_s": sum(int(r["ns"]) for r in block) / 1e9,
                    }
                )
    table(a.out, grouped)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("archived", "structures", "results", "out"):
        p.add_argument(k, type=Path)
    main(p.parse_args())
