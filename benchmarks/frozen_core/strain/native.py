"""Fresh-seed native-search confirmation at fixed pose and affinity selection policies."""

import argparse
import itertools
import statistics
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from pose import GEOM, HERE, physical, read, score_pose, unpack, write


def one(job):
    work, binary, comp, reference, entries = job
    prepared = work / "prepared" / (comp + ".pdbqt")
    scale = physical.evaluate_coordinates(
        binary, prepared, work / "prepared" / (entries[0][2] + ".pdbqt"), None, 90
    )["scale"]
    free8, free32, prep = map(
        float, (reference[k] for k in ["free_8", "free_32", "prepared_mmff"])
    )
    states = []
    for target, assay, rec, seed in entries:
        path = work / "pools" / f"{rec}_{comp}_{seed}.pdbqt"
        poses = physical.read_poses(path)
        physical.fixed_coordinate_inputs(poses, prepared, None)
        for p in poses:
            scores = score_pose(p, reference["smiles"], free8, 0, scale)
            states.append(
                {
                    "target": target,
                    "assay": assay,
                    "receptor": rec,
                    "seed": seed,
                    "compound": comp,
                    "rank": p.native_rank,
                    "native": scores["native"],
                    "unit_8": scores["full_unit"],
                    "unit_32": scores["full_unit"] + scale * (free8 - free32),
                    "prepared_prior_8": scores["native"] + scale * (prep - free8),
                    "prepared_prior_32": scores["native"] + scale * (prep - free32),
                }
            )
    return states


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument(
        "--archive", type=Path, default=HERE.parent / "native/outputs.tar.xz"
    )
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()
    unpack(GEOM / "affinity/inputs.tar.xz", a.work)
    if a.archive:
        unpack(a.archive, a.work)
    ref = {r["compound"]: r for r in read(HERE / "free_references.csv")}
    needed = defaultdict(set)
    for r in read(GEOM / "affinity/cap60000/cells.csv"):
        if r["arm"] == "armijo":
            for seed in (96431, 96457):
                needed[r["compound"]].add(
                    (r["target"], r["assay"], r["receptor"], seed)
                )
    with ProcessPoolExecutor(a.jobs) as pool:
        states = [
            r
            for rows in pool.map(
                one,
                [
                    (a.work, a.binary.resolve(), c, ref[c], sorted(entries))
                    for c, entries in sorted(needed.items())
                ],
            )
            for r in rows
        ]
    write(a.work / "mmff_modes.csv", states)
    preds = defaultdict(list)
    for k, rows in itertools.groupby(
        states,
        key=lambda r: (
            r["target"],
            r["assay"],
            r["receptor"],
            r["compound"],
            r["seed"],
        ),
    ):
        modes = list(rows)
        for arm, policy in itertools.product(
            ["native", "unit_8", "unit_32", "prepared_prior_8", "prepared_prior_32"],
            ["first", "all"],
        ):
            preds[*k[:4], arm, policy].append(
                min(r[arm] for r in modes if policy == "all" or r["rank"] == 1)
            )
    edges = []
    for r in read(GEOM / "affinity/cap60000/pairs.csv"):
        if float(r["fold"]) < 3:
            continue
        k = r["target"], r["assay"], r["receptor"]
        for arm, policy in itertools.product(
            ["native", "unit_8", "unit_32", "prepared_prior_8", "prepared_prior_32"],
            ["first", "all"],
        ):
            weak, strong = (
                preds[*k, r[t], arm, policy] for t in ("weaker", "stronger")
            )
            if len(weak) != 2 or len(strong) != 2:
                raise ValueError("Missing replicate")
            gap = statistics.median(weak) - statistics.median(strong)
            edges.append(
                {
                    "target": k[0],
                    "assay": k[1],
                    "receptor": k[2],
                    "stronger": r["stronger"],
                    "weaker": r["weaker"],
                    "arm": arm,
                    "policy": policy,
                    "gap": gap,
                    "credit": 0.5 if abs(gap) < 1e-10 else float(gap > 0),
                }
            )
    write(a.work / "mmff_pairs.csv", edges)
    groups = defaultdict(list)
    for r in edges:
        groups[r["target"], r["assay"], r["receptor"], r["arm"], r["policy"]].append(
            r["credit"]
        )
    write(
        a.work / "mmff_summary.csv",
        [
            {
                "target": k[0],
                "assay": k[1],
                "receptor": k[2],
                "arm": k[3],
                "policy": k[4],
                "pairs": len(v),
                "credit": sum(v),
            }
            for k, v in sorted(groups.items())
        ],
    )


if __name__ == "__main__":
    main()
