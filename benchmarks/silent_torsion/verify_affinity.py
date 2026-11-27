"""Check raw hashes, first scores, actual counts and tie-aware ranks independently."""

import csv
import hashlib
from itertools import combinations
from pathlib import Path
import re
import tarfile
from fractions import Fraction
import json

HERE = Path(__file__).resolve().parent


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def verify(directory):
    data = read(directory / "cells.csv")
    compounds = {
        r["compound"]: r
        for r in read(HERE.parent / "pair_gradient/galectin_external/manifest.csv")
    }
    scores = {}
    with tarfile.open(directory / "raw.tar.xz") as tar:
        for cell in data:
            arm, rec, lig = [cell[k] for k in ("arm", "receptor", "compound")]
            prefix = f"runs/{rec}_{lig}_{arm}"
            pose = tar.extractfile(prefix + ".pdbqt").read()
            assert hashlib.sha256(pose).hexdigest() == cell["sha256"]
            first = Fraction(
                re.search(rb"REMARK VINA RESULT:\s+([-.\d]+)", pose)[1].decode()
            )
            assert first == Fraction(cell["score"])
            assert (arm, rec, lig) not in scores
            scores[arm, rec, lig] = first
            log = (
                tar.extractfile(prefix + (".counts" if arm == "upstream" else ".log"))
                .read()
                .decode()
            )
            for field in (
                "hunt_evals",
                "second_evals",
                "explicit_evals",
                "hunt_calls",
                "second_calls",
                "explicit_calls",
            ):
                assert int(re.search(r"\b" + field + r"=(\d+)", log)[1]) == int(
                    cell[field]
                )
            cost = json.load(tar.extractfile(prefix + ".cost"))
            assert all(float(cell[k]) == v for k, v in cost.items())
    summary = read(directory / "summary.csv")
    for row in summary:

        def value(lig):
            if row["receptor"] != "mean":
                return scores[row["arm"], row["receptor"], lig]
            return sum(scores[row["arm"], r, lig] for r in ("1KJL", "3T1L")) / 2

        n, credit = 0, 0
        for a, b in combinations(compounds, 2):
            ka, kb = [Fraction(compounds[k]["kd_uM"]) for k in (a, b)]
            if ka == kb or (row["subset"] != "all" and max(ka, kb) < 3 * min(ka, kb)):
                continue
            if row["subset"] == "3fold_matched" and any(
                compounds[a][k] != compounds[b][k] for k in ("heavy", "charge")
            ):
                continue
            pred = (value(b) - value(a)) * (kb - ka)
            n += 1
            credit += (pred > 0) + Fraction(pred == 0, 2)
        assert n == int(row["pairs"]) and credit == Fraction(row["credit"]), row
    print(directory.name, len(data), "docking outputs,", len(summary), "rank summaries")


if __name__ == "__main__":
    for folder in sorted((HERE / "affinity").iterdir()):
        if folder.is_dir():
            verify(folder)
