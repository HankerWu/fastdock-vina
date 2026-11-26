"""Apply the fixed ligand-strain replacements to exact, same-assay Ki or Kd pairs."""

import argparse
import hashlib
import statistics
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from pose import GEOM, METHODS, physical, read, score_pose, unpack, write


def run_ligand(task):
    binary, fix, compound, entries = task
    reference = physical.read_poses(entries[0][-1])[0]
    mol, identity = physical.decode_ligand(reference)
    free, n = physical.gas_reference(mol)
    ligand = fix / (compound + ".pdbqt")
    scale = physical.evaluate_coordinates(
        binary, ligand, fix / (entries[0][1]["receptor"] + ".pdbqt"), None, 90
    )["scale"]
    states, cells = [], []
    for cap, info, path in entries:
        poses = physical.read_poses(path)
        physical.fixed_coordinate_inputs(poses, ligand, None)
        shared = {
            "cap": cap,
            "target": info["target"],
            "assay": info["assay"],
            "receptor": info["receptor"],
            "compound": compound,
            "seed": info["seed"],
        }
        scores = [
            dict(**shared, **score_pose(p, identity, free, n, scale)) for p in poses
        ]
        if min(row["native"] for row in scores) != float(info["top1_score"]):
            raise ValueError("Archived rank one differs")
        cells.append(
            dict(
                **shared,
                pool_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                **{m: min(row[m] for row in scores) for m in METHODS},
            )
        )
        states.extend(scores)
    return states, cells


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()
    source = GEOM / "affinity"
    unpack(source / "inputs.tar.xz", a.work)
    jobs = defaultdict(list)
    for cap in (12000, 60000):
        folder = source / f"cap{cap}"
        local = a.work / str(cap)
        unpack(folder / "outputs.tar.xz", local)
        for r in read(folder / "cells.csv"):
            if r["arm"] == "armijo":
                comp, rec, seed = r["compound"], r["receptor"], r["seed"]
                jobs[comp].append(
                    (cap, r, local / "pools" / f"{rec}_{comp}_{seed}_armijo.pdbqt")
                )
    with ProcessPoolExecutor(a.jobs) as pool:
        results = list(
            pool.map(
                run_ligand,
                [
                    (a.binary.resolve(), a.work / "prepared", comp, v)
                    for comp, v in sorted(jobs.items())
                ],
            )
        )
    states = [r for rows, _ in results for r in rows]
    cells = [r for _, rows in results for r in rows]
    write(a.work / "modes.csv", states)
    write(a.work / "cells.csv", cells)
    predicted = defaultdict(list)
    for r in cells:
        for method in METHODS:
            predicted[
                r["cap"], r["target"], r["assay"], r["receptor"], r["compound"], method
            ].append(r[method])
    pairs = []
    for cap in (12000, 60000):
        for r in read(source / f"cap{cap}" / "pairs.csv"):
            if float(r["fold"]) < 3:
                continue
            group = cap, r["target"], r["assay"], r["receptor"]
            for method in METHODS:
                scores = [
                    predicted[*group, r[k], method] for k in ("weaker", "stronger")
                ]
                if not all(len(v) == 2 for v in scores):
                    raise ValueError("Expected two seeds per ligand")
                delta = statistics.median(scores[0]) - statistics.median(scores[1])
                acc = 0.5 if abs(delta) < 1e-10 else float(delta > 0)
                if method == "native" and (
                    abs(delta - float(r["armijo_delta"])) > 1e-7
                    or acc != float(r["armijo_correct"])
                ):
                    raise ValueError("Mismatch with original assay concordance")
                pairs.append(
                    {
                        "cap": cap,
                        "target": r["target"],
                        "assay": r["assay"],
                        "receptor": r["receptor"],
                        "stronger": r["stronger"],
                        "weaker": r["weaker"],
                        "fold": r["fold"],
                        "method": method,
                        "weak_minus_strong": delta,
                        "credit": acc,
                    }
                )
    write(a.work / "pairs.csv", pairs)
    groups = defaultdict(list)
    for r in pairs:
        groups[r["cap"], r["target"], r["assay"], r["receptor"], r["method"]].append(
            r["credit"]
        )
    write(
        a.work / "summary.csv",
        [
            {
                "cap": k[0],
                "target": k[1],
                "assay": k[2],
                "receptor": k[3],
                "method": k[4],
                "pairs": len(v),
                "credit": sum(v),
            }
            for k, v in sorted(groups.items())
        ],
    )


if __name__ == "__main__":
    main()
