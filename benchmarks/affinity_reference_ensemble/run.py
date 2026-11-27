"""Apply the fixed protocol to archived native top poses from both Kd studies."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import tarfile
import time

import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
CONFIG = json.loads((HERE / "protocol.json").read_text())
STATES = {
    "1GZW": "Gal1",
    "4Q26": "Gal1",
    "1KJL": "Gal3",
    "3T1L": "Gal3",
    "2YWP": "CHEK1",
    "2E9P": "CHEK1",
    "2J4Z": "AURKA",
    "2X6E": "AURKA",
}


def direction(u, v):
    z = 2 * u - 1
    radius = np.sqrt(1 - z * z)
    phi = 2 * np.pi * v
    return np.column_stack([radius * np.cos(phi), radius * np.sin(phi), z])


def nodes(path, count=None, scramble=40721):
    n = count or max(CONFIG["convergence_points"] + [CONFIG["primary_points"]])
    u = qmc.Sobol(d=7, scramble=True, seed=scramble).random_base2(int(np.log2(n)))
    shift = (
        direction(u[:, 0], u[:, 1])
        * np.cbrt(u[:, 2, None])
        * CONFIG["translation_radius_A"]
    )
    cap = np.deg2rad(CONFIG["rotation_half_angle_deg"])
    target = u[:, 5] * (cap - np.sin(cap))
    lo, hi = np.zeros(n), np.full(n, cap)
    for _ in range(48):
        mid = 0.5 * (lo + hi)
        smaller = mid - np.sin(mid) < target
        lo, hi = np.where(smaller, mid, lo), np.where(smaller, hi, mid)
    rot = Rotation.from_rotvec(
        direction(u[:, 3], u[:, 4]) * ((lo + hi) / 2)[:, None]
    ).as_matrix()
    with path.open("w") as f:
        np.savetxt(
            f, np.column_stack([shift, rot.reshape(-1, 9), u[:, 6]]), fmt="%.16g"
        )


def extract(directory, archive):
    expected = (archive.parent / "evidence.sha256").read_text().split()[0]
    if hashlib.sha256(archive.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Wrong raw inputs: {archive}")
    with tarfile.open(archive) as tf:
        for entry in tf:
            p = Path(entry.name)
            if (
                len(p.parts) == 2
                and p.suffix == ".pdbqt"
                and (
                    p.parent.name == "prepared"
                    or (p.parent.name == "pools" and p.name.endswith("_native.pdbqt"))
                )
            ):
                destination = directory / p
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(tf.extractfile(entry).read())


def labels(repo):
    result = {}
    for row in csv.DictReader(
        (repo / "benchmarks/galectin_affinity/table1.csv").open()
    ):
        for target in ("Gal1", "Gal3"):
            result[target, row["compound"]] = float(row[target + "_Kd_uM"]) * 1000
    for row in csv.DictReader(
        (repo / "benchmarks/torsion_kinase_affinity/panel.csv").open()
    ):
        if row["compound"] != "Ki-20227":
            result[row["target"], row["cid"]] = float(row["kd_nm"])
    return result


def run_one(arguments):
    binary, points, inp, kd, *bins = arguments
    state, ligand, seed, _ = inp.stem.split("_")
    prepared = inp.parent.parent / "prepared"
    proc = subprocess.run(
        [
            binary,
            prepared / (state + ".pdbqt"),
            prepared / (ligand + ".pdbqt"),
            inp,
            points,
            str(CONFIG["temperature_K"]),
            *map(str, bins),
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    lines = inp.read_text().split("ENDMDL", 1)[0].splitlines()

    def number(prefix):
        return float(
            next(s for s in lines if s.startswith(prefix))[len(prefix) :].split()[0]
        )

    named = dict(
        target=STATES[state],
        state=state,
        ligand=ligand,
        seed=seed,
        kd_nM=kd,
        native=number("REMARK VINA RESULT:"),
        printed_inter=number("REMARK INTER:"),
    )
    return [dict(**named, **row) for row in csv.DictReader(io.StringIO(proc.stdout))]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("repo", type=Path)
    p.add_argument("integrator", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--nodes", type=int, default=max(CONFIG["convergence_points"]))
    p.add_argument("--scramble", type=int, default=40721)
    p.add_argument("--bins", type=int, default=16)
    a = p.parse_args()
    if a.nodes < 256 or a.nodes & (a.nodes - 1):
        p.error("--nodes must be a power of two >=256")
    a.work.mkdir(parents=True, exist_ok=True)
    points = a.work / "nodes.txt"
    nodes(points, a.nodes, a.scramble)
    target_kd = labels(a.repo)
    todo = []
    for source in ("galectin_affinity", "torsion_kinase_affinity"):
        directory = a.work / source
        extract(directory, a.repo / "benchmarks" / source / "evidence.tar.xz")
        for pool in sorted((directory / "pools").glob("*_native.pdbqt")):
            state, ligand, _, _ = pool.stem.split("_")
            key = (STATES[state], ligand)
            if key in target_kd:
                todo.append(
                    (a.integrator.resolve(), points, pool, target_kd[key], a.bins)
                )
    if len(todo) != 8 * len(target_kd):
        raise ValueError("Missing native receptor/seed blocks")
    began = time.monotonic()
    before = os.times()
    with (
        ThreadPoolExecutor(max_workers=a.workers) as workers,
        (a.work / "cells.csv").open("w") as f,
    ):
        writer = None
        for result in workers.map(run_one, todo):
            if writer is None:
                writer = csv.DictWriter(f, fieldnames=result[0].keys())
                writer.writeheader()
            writer.writerows(result)
            f.flush()
    after = os.times()
    (a.work / "timing.json").write_text(
        json.dumps(
            dict(
                workers=a.workers,
                cell_count=len(todo),
                wall_s=time.monotonic() - began,
                child_cpu_s=after.children_user
                + after.children_system
                - before.children_user
                - before.children_system,
                scramble=a.scramble,
                hinge_bins=a.bins,
                joint_nodes=a.nodes,
                rigid_control_nodes=a.nodes,
                hinge_control_nodes=a.bins,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
