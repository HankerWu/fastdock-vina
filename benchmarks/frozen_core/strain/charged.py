"""Post-hoc dielectric sensitivity for the two +1 BACE1 Ki ligand series."""

import argparse
import statistics
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from pose import ROOT, physical, read, score_pose, write

sys.path.insert(0, str(ROOT / "benchmarks/strain_transfer"))
from evaluate import screened_properties, screened_reference


def one(job):
    work, comp, cond = job
    ligand = work / "prepared" / (comp + ".pdbqt")
    initial, identity = physical.decode_ligand(
        physical.Pose(0, 0, 0, ligand.read_text(), None, None)
    )
    free, nfree = screened_reference(initial)
    rows = []
    for key, scale in cond.items():
        cap, _, rec, seed = key
        path = work / cap / "pools" / f"{rec}_{comp}_{seed}_armijo.pdbqt"
        poses = physical.read_poses(path)
        physical.fixed_coordinate_inputs(poses, ligand, None)
        values = [
            score_pose(p, identity, free, nfree, scale, screened_properties)
            for p in poses
        ]
        best = min(values, key=lambda m: (m["full_unit"], m["rank"]))
        rows.append(
            {
                "cap": cap,
                "assay": key[1],
                "receptor": rec,
                "compound": comp,
                "seed": seed,
                "screened_rank": best["rank"],
                "screened_bound": best["fixed_mmff"],
                "screened_free": free,
                "screened_score": best["full_unit"],
            }
        )
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("affinity_work", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()
    cond = defaultdict(dict)
    for r in read(a.affinity_work / "modes.csv"):
        if r["target"] == "BACE1":
            cond[r["compound"]][r["cap"], r["assay"], r["receptor"], r["seed"]] = float(
                r["flexibility_scale"]
            )
    with ProcessPoolExecutor(a.jobs) as pool:
        rows = [
            r
            for rs in pool.map(
                one, [(a.affinity_work, comp, v) for comp, v in sorted(cond.items())]
            )
            for r in rs
        ]
    write(a.affinity_work / "screened_cells.csv", rows)
    by = defaultdict(list)
    for r in rows:
        by[r["cap"], r["assay"], r["receptor"], r["compound"]].append(
            r["screened_score"]
        )
    scores = {k: statistics.median(v) for k, v in by.items()}
    summary = defaultdict(list)
    for r in read(a.affinity_work / "pairs.csv"):
        if r["target"] == "BACE1" and r["method"] == "native":
            block = r["cap"], r["assay"], r["receptor"]
            gap = scores[*block, r["weaker"]] - scores[*block, r["stronger"]]
            credit = 0.5 if abs(gap) < 1e-10 else float(gap > 0)
            summary[block].append(credit)
    write(
        a.affinity_work / "screened_summary.csv",
        [
            {
                "cap": k[0],
                "assay": k[1],
                "receptor": k[2],
                "pairs": len(v),
                "credit": sum(v),
            }
            for k, v in sorted(summary.items())
        ],
    )


if __name__ == "__main__":
    main()
