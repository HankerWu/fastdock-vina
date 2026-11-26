"""Hindsight conditional pose selection on fixed transferred pools; scores are within-ligand only."""

import argparse
import csv
import re
from pathlib import Path

import numpy as np
from transfer_holdout import write_csv


def evaluate(runs):
    with (runs / "modes.csv").open(newline="") as stream:
        original = {
            (r["ligand"], r["receptor"], r["seed"], int(r["model"])): float(r["rmsd"])
            for r in csv.DictReader(stream)
            if r["arm"] == "transfer"
        }
    modes, winners = [], []
    for path in sorted(runs.glob("*_transfer.pdbqt")):
        rec, lig, seed, _ = path.stem.split("_")
        with (runs / f"{rec}_{lig}.anchor").open() as stream:
            anchors = [line.split() for line in stream if line.strip()]
        if not anchors or any(len(a) != 4 for a in anchors):
            raise ValueError(f"Point-only transfer prior expected: {rec}/{lig}")
        serials = [int(a[0]) for a in anchors]
        targets = np.array([[float(x) for x in a[1:]] for a in anchors])
        cell = []
        for model, block in enumerate(
            re.findall(
                r"^MODEL .*?^ENDMDL$", path.read_text(), re.MULTILINE | re.DOTALL
            ),
            1,
        ):
            xyz = {
                int(l[6:11]): [float(l[j : j + 8]) for j in (30, 38, 46)]
                for l in block.splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
            }
            offsets = np.array([xyz[serial] for serial in serials]) - targets
            squared = np.sum(offsets**2, axis=1)
            penalty = float(np.mean(4 * squared / (np.sqrt(1 + squared) + 1)))
            raw = 0.0
            for term in ("INTER", "INTRA"):
                raw += float(
                    re.search(r"^REMARK " + term + r":\s+(\S+)", block, re.MULTILINE)[1]
                )
            r = {
                "receptor": rec,
                "ligand": lig,
                "seed": seed,
                "model": model,
                "rmsd": original[(lig, rec, seed, model)],
                "raw": raw,
                "penalty": penalty,
                "composite": raw + penalty,
            }
            cell.append(r)
        modes.extend(cell)
        for ranker in ("released", "composite", "penalty"):
            ordered = (
                cell
                if ranker == "released"
                else sorted(cell, key=lambda x: (x[ranker], x["model"]))
            )
            top = ordered[0]
            winners.append(
                dict(
                    **top,
                    ranker=ranker,
                    first_near=next(
                        (r for r, m in enumerate(ordered, 1) if m["rmsd"] <= 2), None
                    ),
                    oracle=min(r["rmsd"] for r in cell),
                )
            )
    return modes, winners


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("runs", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    m, w = evaluate(a.runs)
    write_csv(a.output / "conditional_candidates.csv", m)
    write_csv(a.output / "conditional_top.csv", w)
