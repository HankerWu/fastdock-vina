"""Compare matched BRD4 native-scored poses and Kd order across three search arms."""

import argparse
import csv
import itertools
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median
from tempfile import TemporaryDirectory

import numpy as np
from experiment import ARMS, LIGANDS, RECEPTORS, ROOT, SEEDS, write

sys.path.insert(0, str(ROOT / "benchmarks/alkyl_holdout"))
from evaluate import read_poses

sys.path.insert(0, str(ROOT / "benchmarks/cutoff_switch"))
from verify_cutoff_poses import verify


def summarize(cells, pairs):
    summary = []
    for lig, arm in itertools.product(("all", *LIGANDS), ARMS):
        block = [
            x for x in cells if x["arm"] == arm and (lig == "all" or x["ligand"] == lig)
        ]
        p = [x for x in pairs if x["arm"] == arm]
        summary.append(
            {
                "ligand": lig,
                "arm": arm,
                "n": len(block),
                "top2A": sum(x["top_rmsd"] <= 2 for x in block),
                "pool2A": sum(x["pool_rmsd"] <= 2 for x in block),
                "median_top_rmsd": median(x["top_rmsd"] for x in block),
                "pairs": len(p) if lig == "all" else "",
                "pair_credit": sum(x["credit"] for x in p) if lig == "all" else "",
            }
        )
    return summary


def evaluate(data, work, controls, out):
    out.mkdir(exist_ok=True, parents=True)
    meta = json.loads((data / "prepared/manifest.json").read_text())
    sites = defaultdict(list)
    for row in csv.DictReader((work / "sites.csv").open()):
        sites[row["receptor"], row["ligand"]].append(
            (int(row["serial"]), np.array([float(row[k]) for k in "xyz"]))
        )
    modes = []
    cells = []
    pairs = []
    geometry = {}
    with TemporaryDirectory() as tmp:
        audit = Path(tmp)
        (audit / "raw").symlink_to(data.resolve() / "raw")
        (audit / "runs").mkdir()
        comparisons = []
        for rec, seed, arm in itertools.product(RECEPTORS, SEEDS, ARMS):
            scores = {}
            for lig in LIGANDS:
                stem = (
                    work if arm == "contact" else controls
                ) / f"{rec}_{lig}_{seed}_{arm}"
                log = stem.with_suffix(".log").read_text()
                evs = [
                    int(x)
                    for x in re.findall(r"(?:hunt_evals|second_evals)=(\d+)", log)
                ]
                if len(evs) != 8 or not all(
                    a + b >= 50000 for a, b in zip(evs[::2], evs[1::2], strict=True)
                ):
                    raise ValueError(f"Incomplete search {stem}")
                if (
                    arm == "contact"
                    and f"sites=1 assignments={len(sites[rec, lig])}" not in log
                ):
                    raise ValueError("Incorrect alternative count")
                src = stem.with_suffix(".pdbqt").resolve()
                (audit / "runs" / f"{lig}_{rec}_{seed}_{arm}.pdbqt").symlink_to(src)
                serials = meta[lig]["heavy_serials"]
                refs = np.asarray(meta[lig]["alignments"][rec]["reference_by_match"])
                points = sites[rec, lig]
                geo = min(
                    np.linalg.norm(refs[:, serials.index(s), :] - target, axis=1).min()
                    for s, target in points
                )
                geometry[rec, lig] = {
                    "receptor": rec,
                    "ligand": lig,
                    "min_crystal_to_site_A": geo,
                    "assignment_count": len(points),
                }
                base = {"receptor": rec, "ligand": lig, "seed": seed, "arm": arm}
                pool = []
                for rank, score, coords in read_poses(src):
                    xyz = np.array([coords[s] for s in serials])
                    rms = np.sqrt(np.square(refs - xyz).sum(axis=2).mean(axis=1)).min()
                    err = min(
                        np.linalg.norm(np.array(coords[s]) - t) for s, t in points
                    )
                    r = {
                        **base,
                        "rank": rank,
                        "score": score,
                        "rmsd": rms,
                        "nearest_contact_site_A": err,
                    }
                    modes.append(r)
                    pool.append(r)
                if not 1 <= len(pool) <= 9 or pool[0]["rank"] != 1:
                    raise ValueError(f"Unexpected returned models: {src}")
                cells.append(
                    {
                        **base,
                        "top_rmsd": pool[0]["rmsd"],
                        "pool_rmsd": min(x["rmsd"] for x in pool),
                        "top_score": pool[0]["score"],
                        "grid_evals": sum(evs),
                        "n_models": len(pool),
                    }
                )
                scores[lig] = pool[0]["score"]
                comparisons.append(
                    {**base, "pool_rmsds": ",".join(str(x["rmsd"]) for x in pool)}
                )
            ordered = sorted(LIGANDS, key=lambda l: meta[l]["kd_uM"])
            for strong, weak in itertools.combinations(ordered, 2):
                ratio = meta[weak]["kd_uM"] / meta[strong]["kd_uM"]
                if ratio >= 3:
                    gap = scores[weak] - scores[strong]
                    pairs.append(
                        {
                            "receptor": rec,
                            "seed": seed,
                            "arm": arm,
                            "strong": strong,
                            "weak": weak,
                            "kd_ratio": ratio,
                            "score_margin": gap,
                            "credit": 1 if gap > 0 else 0.5 if gap == 0 else 0,
                        }
                    )
        write(audit / "cells.csv", comparisons)
        verify(audit)
    summary = summarize(cells, pairs)
    for name, rows in (
        ("modes", modes),
        ("cells", cells),
        ("pairs", pairs),
        ("summary", summary),
        ("geometry", list(geometry.values())),
    ):
        write(out / f"{name}.csv", rows)
    print((out / "summary.csv").read_text())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("controls", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    evaluate(a.data, a.work, a.controls, a.out)
