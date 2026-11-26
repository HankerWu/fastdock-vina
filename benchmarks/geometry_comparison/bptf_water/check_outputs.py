"""Recompute water-arm ranks, scores and RMSDs from archived pose coordinates."""

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pose_run import align_with_numbering, load_meta


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fixtures", type=Path)
    p.add_argument("work", type=Path)
    args = p.parse_args()
    meta = {pdb: load_meta(args.fixtures, pdb) for pdb in ("7F5D", "7RWO")}
    for row in json.loads((args.work / "external_manifest.json").read_text()):
        if row["status"] == "ready":
            meta[row["pdb"]] = load_meta(args.work / "external", row["pdb"])
    for prefix, subdir in (("", "poses"), ("confirm_", "confirm")):
        cells = list(csv.DictReader((args.work / (prefix + "cells.csv")).open()))
        mode_rows = defaultdict(list)
        for row in csv.DictReader((args.work / (prefix + "modes.csv")).open()):
            mode_rows[row["query"], row["receiver"], row["seed"], row["arm"]].append(
                row
            )
        count = 0
        for row in cells:
            q, r, seed, arm = (row[key] for key in ("query", "receiver", "seed", "arm"))
            path = args.work / arm / subdir / f"{q}-{r}-{seed}-{arm}.pdbqt"
            assert hashlib.sha256(path.read_bytes()).hexdigest() == row["pdbqt_sha256"]
            ref = np.asarray(meta[q]["mapped"])
            if q != r:
                rot, tr, *_ = align_with_numbering(meta[q], meta[r])
                ref = ref @ rot + tr
            values = []
            for block in path.read_text().split("ENDMDL"):
                xyz = {
                    int(l[6:11]): [float(l[i : i + 8]) for i in (30, 38, 46)]
                    for l in block.splitlines()
                    if l.startswith(("ATOM  ", "HETATM"))
                }
                if not xyz:
                    continue
                rms = min(
                    math.sqrt(
                        sum(
                            math.dist(xyz[i], pt) ** 2
                            for i, pt in zip(meta[q]["serial"], m)
                        )
                        / len(m)
                    )
                    for m in ref
                )
                expected = mode_rows[q, r, seed, arm][len(values)]
                assert abs(rms - float(expected["rmsd_A"])) < 1e-8
                score = next(
                    float(l.split()[3])
                    for l in block.splitlines()
                    if l.startswith("REMARK VINA RESULT:")
                )
                assert score == float(expected["score"])
                values.append(rms)
            assert len(values) == int(row["modes"])
            assert abs(min(values) - float(row["oracle_A"])) < 1e-8
            assert abs(values[0] - float(row["top1_A"])) < 1e-8
            count += len(values)
        print(prefix or "pilot", len(cells), "cells", count, "modes verified")


if __name__ == "__main__":
    main()
