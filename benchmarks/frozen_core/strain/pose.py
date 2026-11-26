"""Frozen-pool, cross-target transfer test of full MMFF94s intraligand strain."""

import argparse
import csv
import hashlib
import json
import sys
import tarfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
GEOM = ROOT / "benchmarks/geometry_comparison"
sys.path.insert(0, str(ROOT / "tools"))
import score_physical_pool as physical

sys.path.insert(0, str(GEOM / "upstream"))
import replay as upstream

METHODS = ("native", "inter_only", "historical_01", "full_unit")


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def unpack(path, dest):
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path) as archive:
        archive.extractall(dest, filter="data")


def collect(work):
    inputs = work / "inputs"
    unpack(GEOM / "pose/inputs.tar.xz", inputs / "cognate")
    unpack(
        ROOT / "benchmarks/cubic_grid/human_panel/prepared.tar.xz",
        inputs / "reciprocal",
    )
    jobs = defaultdict(list)
    for budget in ["low", "high"]:
        dest = work / budget
        unpack(GEOM / ("pose/" + budget + "_outputs.tar.xz"), dest)
        modes = defaultdict(list)
        for row in read(GEOM / ("pose/" + budget + "_modes.csv")):
            if row["arm"] == "armijo":
                modes[row["query"], row["receiver"], row["seed"]].append(
                    float(row["rmsd_A"])
                )
        for row in read(GEOM / ("pose/" + budget + "_cells.csv")):
            if row["arm"] != "armijo":
                continue
            q, r, seed = row["query"], row["receiver"], row["seed"]
            path = dest / "poses" / f"{q}-{r}-{seed}-armijo.pdbqt"
            jobs[q].append(
                (
                    row["target"],
                    r,
                    seed,
                    budget,
                    inputs / "cognate",
                    path,
                    row["pdbqt_sha256"],
                    modes[q, r, seed],
                )
            )

    for protocol in ["cognate", "reciprocal"]:
        fixtures = inputs / protocol
        framed = {
            (q["pdb"], r["pdb"]): (q, ref)
            for q, r, ref in upstream.experiments(protocol, fixtures)
        }
        for budget, arm in [("repeat", "armijo"), ("matched", "native")]:
            src = GEOM / "upstream" / protocol
            if budget == "matched":
                src = src / "matched"
            name = budget + "_" + protocol
            dest = work / name
            unpack(src / "raw.tar.xz", dest)
            for row in read(src / "cells.csv"):
                if row["arm"] != arm:
                    continue
                q, r, seed = row["query"], row["receiver"], row["seed"]
                path = dest / f"{q}-{r}-{seed}-{arm}.pdbqt"
                qm, expected = framed[q, r]
                rmsds = upstream.reciprocal.measure(path, qm["serial"], expected)
                if abs(rmsds[0] - float(row["first"])) > 1e-6:
                    raise ValueError("Archived native pose differs")
                jobs[q].append(
                    (
                        row["target"],
                        r,
                        seed,
                        name,
                        fixtures,
                        path,
                        row["pool_sha256"],
                        rmsds,
                    )
                )
    return jobs


def score_pose(
    p, identity, free_E, converged, scale, properties=physical.mmff_properties
):
    mol, current = physical.decode_ligand(p)
    if current != identity:
        raise ValueError("Ligand graph or stereochemistry changed")
    conf = mol.GetConformer()
    idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    pos = conf.GetPositions()[idx]
    energy = physical.minimized_energy(mol, properties(mol), conf.GetId(), idx)
    if not np.array_equal(pos, conf.GetPositions()[idx]):
        raise ValueError("Heavy-atom coordinates moved")
    native = p.native_score
    inter = native - scale * (p.native_intra - p.native_unbound)
    delta = scale * (energy - free_E)
    return {
        "rank": p.model,
        "vina_intra": p.native_intra,
        "vina_unbound": p.native_unbound,
        "flexibility_scale": scale,
        "fixed_mmff": energy,
        "free_mmff": free_E,
        "nfree": converged,
        "native": native,
        "inter_only": inter,
        "historical_01": inter + 0.1 * delta,
        "full_unit": physical.mmff_affinity_score(p, scale, energy, free_E),
    }


def process(task):
    binary, q, entries = task
    first = physical.read_poses(entries[0][5])[0]
    refmol, identity = physical.decode_ligand(first)
    free_E, converged = physical.gas_reference(refmol)
    fixes = entries[0][4]
    scale = physical.evaluate_coordinates(
        binary, fixes / q / "ligand.pdbqt", fixes / q / "receptor.pdbqt", None, 90
    )["scale"]
    modes, cells = [], []
    for target, receiver, seed, panel, fix, path, digest, rmsds in entries:
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError("Frozen pool checksum mismatch " + str(path))
        poses = physical.read_poses(path)
        physical.fixed_coordinate_inputs(poses, fix / q / "ligand.pdbqt", None)
        if len({p.native_unbound for p in poses}) != 1:
            raise ValueError("Expected a constant unbound reference")
        shared = {
            "panel": panel,
            "target": target,
            "query": q,
            "receiver": receiver,
            "seed": seed,
        }
        pool = []
        for p, rmsd in zip(poses, rmsds, strict=True):
            pool.append(
                dict(
                    **shared,
                    rmsd_A=rmsd,
                    **score_pose(p, identity, free_E, converged, scale),
                )
            )
        modes.extend(pool)
        for method in METHODS:
            ordered = sorted(pool, key=lambda row: (row[method], row["rank"]))
            cells.append(
                dict(
                    **shared,
                    method=method,
                    native_rank_of_top=ordered[0]["rank"],
                    top1_A=ordered[0]["rmsd_A"],
                    top5_A=min(m["rmsd_A"] for m in ordered[:5]),
                    oracle_A=min(m["rmsd_A"] for m in pool),
                    modes=len(pool),
                    pool_sha256=digest,
                )
            )
    return modes, cells


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    args = p.parse_args()
    work = args.work.resolve()
    jobs = collect(work)
    with ProcessPoolExecutor(args.jobs) as pool:
        results = list(
            pool.map(
                process,
                [
                    (args.binary.resolve(), q, entries)
                    for q, entries in sorted(jobs.items())
                ],
            )
        )
    rows = [v for modes, _ in results for v in modes]
    cells = [v for _, out in results for v in out]
    write(work / "modes.csv", rows)
    write(work / "cells.csv", cells)
    summary = []
    for label in sorted({v["panel"] for v in cells}):
        for m in METHODS:
            subset = [v for v in cells if v["panel"] == label and v["method"] == m]
            summary.append(
                {
                    "panel": label,
                    "method": m,
                    "pools": len(subset),
                    "top1": sum(v["top1_A"] <= 2 for v in subset),
                    "top5": sum(v["top5_A"] <= 2 for v in subset),
                    "oracle": sum(v["oracle_A"] <= 2 for v in subset),
                }
            )
    write(work / "summary.csv", summary)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
