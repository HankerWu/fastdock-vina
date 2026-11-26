"""Fixed-pose rotor integration on the immutable native Late pools."""

import argparse
import csv
import io
import itertools
import math
import statistics
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmarks/strain_transfer"))
from evaluate import SHARED, load, panels, write

# T4L's 29 C calorimetry; other targets use a declared reference temperature.
TEMPERATURE = {"t4l": 302.15, "brd4": 298.15, "hcaii": 298.15}


def run_cell(job, binary, evidence, angles):
    (target, lig), cell = job
    rec = cell["receptor"]
    prep = evidence / target / "prepared"
    args = [
        binary,
        prep / f"{rec}_receptor.pdbqt",
        prep / f"{lig}_ligand.pdbqt",
        cell["path"],
        str(angles),
        str(TEMPERATURE[target]),
    ]
    result = subprocess.run(args, capture_output=True, text=True, check=True)
    scores = list(csv.DictReader(io.StringIO(result.stdout)))
    rmsds = list(map(float, cell["pool_rmsds"].split(",")))
    output = []
    for r, rms in zip(scores, rmsds, strict=True):
        output.append(
            dict(
                target=target,
                ligand=lig,
                receptor=rec,
                seed=int(cell["seed"]),
                rmsd=rms,
                measurement=float(cell.get("measurement", cell.get("kd_uM"))),
                model=int(r["model"]),
                rotors=int(r["rotors"]),
                **{
                    k: float(r[k]) for k in ("inter", "intra", "work", "half_grid_work")
                },
            )
        )
    return output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("evidence", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--angles", type=int, default=72)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    labels = load(ROOT / "benchmarks/strain_transfer/results/poses.csv")
    keyed = {
        (r["target"], r["ligand"], r["receptor"], int(r["seed"]), int(r["model"])): r
        for r in labels
    }
    poses = []
    with ThreadPoolExecutor(a.jobs) as pool:
        for block in pool.map(
            lambda c: run_cell(c, a.binary.resolve(), a.evidence, a.angles),
            panels(a.evidence),
        ):
            for r in block:
                k = (r["target"], r["ligand"], r["receptor"], r["seed"], r["model"])
                old = keyed[k]
                r["native"] = float(old["native"])
                scale = float(old["scale"])
                old_inter = (
                    r["native"] / scale
                    - float(old["native_intra"])
                    + float(old["native_unbound"])
                )
                if (
                    abs(r["inter"] - old_inter) > 0.10
                    or abs(r["intra"] - float(old["native_intra"])) > 0.10
                ):
                    raise ValueError(
                        f"Reconstructed energy differs: {k} {r['inter']} {old_inter}"
                    )
                r["no_divisor"] = r["inter"] + r["intra"] - float(old["native_unbound"])
                r["rotor"] = r["no_divisor"] + r["work"]
                r["uniform_three"] = r["no_divisor"] + 0.00198720425864 * TEMPERATURE[
                    r["target"]
                ] * r["rotors"] * math.log(3)
                poses.append(r)
    write(a.output / "rotor_poses.csv", poses)
    cells = []
    by_pool = lambda r: (r["target"], r["receptor"], r["seed"], r["ligand"])
    for (target, rec, seed, lig), group in itertools.groupby(
        sorted(poses, key=by_pool), by_pool
    ):
        group = list(group)
        for method in (
            "native",
            "no_divisor",
            "uniform_three",
            "rotor",
            "rotor_reselect",
        ):
            key = "rotor" if method == "rotor_reselect" else method
            top = (
                min(group, key=lambda x: (x[key], x["model"]))
                if method == "rotor_reselect"
                else group[0]
            )
            assert group[0]["model"] == 1
            cells.append(
                {
                    "target": target,
                    "receptor": rec,
                    "seed": seed,
                    "ligand": lig,
                    "arm": method,
                    "model": top["model"],
                    "rmsd": top["rmsd"],
                    "measurement": top["measurement"],
                    "score": top[key],
                }
            )
    write(a.output / "rotor_cells.csv", cells)
    blocks = defaultdict(list)
    for r in cells:
        if r["receptor"] in SHARED[r["target"]]:
            blocks[(r["target"], r["receptor"], r["seed"], r["arm"])].append(r)
    pairs = []
    for (target, rec, seed, arm), block in sorted(blocks.items()):
        for s, w in itertools.combinations(
            sorted(block, key=lambda x: x["measurement"]), 2
        ):
            if w["measurement"] < 3 * s["measurement"]:
                continue
            gap = w["score"] - s["score"]
            pairs.append(
                {
                    "target": target,
                    "receptor": rec,
                    "seed": seed,
                    "arm": arm,
                    "strong": s["ligand"],
                    "weak": w["ligand"],
                    "gap": gap,
                    "correct": int(gap > 0) + 0.5 * (gap == 0),
                }
            )
    write(a.output / "rotor_pairs.csv", pairs)
    summary = []
    for target, arm in sorted({(p["target"], p["arm"]) for p in pairs}):
        obs = [p for p in pairs if (p["target"], p["arm"]) == (target, arm)]
        unique = defaultdict(list)
        for p in obs:
            unique[(p["strong"], p["weak"])].append(p["gap"])
        geom = [c for c in cells if (c["target"], c["arm"]) == (target, arm)]
        gaps = [statistics.median(v) for v in unique.values()]
        summary.append(
            {
                "target": target,
                "arm": arm,
                "pose_cells": len(geom),
                "top_hits": sum(r["rmsd"] <= 2 for r in geom),
                "credit": sum(p["correct"] for p in obs),
                "total": len(obs),
                "unique_credit": sum(int(g > 0) + 0.5 * (g == 0) for g in gaps),
                "unique_pairs": len(gaps),
            }
        )
    write(a.output / "rotor_summary.csv", summary)


if __name__ == "__main__":
    main()
