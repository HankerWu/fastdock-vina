"""New-target transfer. Rebuild the frozen S6K1 subset, dock it and marginalize each Pro top pose."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

from run import extract, nodes, run_one, STATES

STATES["3A60"] = "S6K1"


def write(path, rows):
    with path.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def prepare(args):
    import xlrd

    raw = args.work / "raw"
    sources = json.loads(
        (args.repo / "benchmarks/kinase_state_affinity/source/sources.json").read_text()
    )
    for source in ("table1.xls", "table4.xls"):
        if (
            hashlib.sha256((raw / source).read_bytes()).hexdigest()
            != sources[source]["sha256"]
        ):
            raise ValueError(f"Unexpected original Davis file: {source}")
    if (
        hashlib.sha256((raw / "3A60.pdb").read_bytes()).hexdigest()
        != "0ace4ecc3ef48435bd8e18ee30639cf11cfec53e306176c4ad519614a382cb14"
    ):
        raise ValueError("Unexpected 3A60 coordinate revision")
    tables = [
        xlrd.open_workbook(str(raw / name)).sheet_by_index(0)
        for name in ("table1.xls", "table4.xls")
    ]
    meta, target = [
        next(row for row in sheet.get_rows() if row[2].value == "S6K1")
        for sheet in tables
    ]
    if (
        meta[3].value != "NO"
        or meta[0].value != target[0].value
        or meta[1].value != "RPS6KB1"
    ):
        raise ValueError("Assay identity mismatch")
    known = {
        r["compound"]: r["cid"]
        for r in csv.DictReader(
            (args.repo / "benchmarks/torsion_kinase_affinity/panel.csv").open()
        )
    }
    rows = []
    for col, name in enumerate(tables[1].row_values(0)):
        if name not in known:
            continue
        value = target[col].value
        reason = (
            "no exact Kd"
            if not value
            else "cocrystal donor"
            if name == "Staurosporine"
            else "racemate"
            if name == "Ki-20227"
            else "included"
        )
        rows.append(dict(compound=name, cid=known[name], kd_nM=value, decision=reason))
    write(args.work / "panel.csv", rows)
    extract(
        args.work / "inputs",
        args.repo / "benchmarks/torsion_kinase_affinity/evidence.tar.xz",
    )
    sys.path.insert(0, str(args.repo / "benchmarks/torsion_kinase_affinity"))
    from prepare_receptors import prepare as receptor

    receptor(args.work, "3A60", "STU")
    for r in rows:
        if r["decision"] == "included":
            shutil.copyfile(
                args.work / "inputs/prepared" / (r["cid"] + ".pdbqt"),
                args.work / "prepared" / (r["cid"] + ".pdbqt"),
            )


def dock(args):
    rows = [
        r
        for r in csv.DictReader((args.work / "panel.csv").open())
        if r["decision"] == "included"
    ]
    center = json.loads((args.work / "prepared/3A60.json").read_text())["center"]
    ref = args.work / "prepared/3A60.pdbqt"
    runs = [
        (r, seed, arm)
        for r in rows
        for seed in (271828, 314159)
        for arm in ("pro", "pro_double", "upstream")
    ]
    pools = args.work / "pools"
    pools.mkdir(exist_ok=True)
    points = args.work / "nodes.txt"
    nodes(points)

    def cell(job):
        r, seed, arm = job
        binary = args.native if arm == "upstream" else args.pro
        pool = pools / f"3A60_{r['cid']}_{seed}_{arm}.pdbqt"
        cmd = [
            str(binary),
            "--receptor",
            str(ref),
            "--ligand",
            str(ref.parent / (r["cid"] + ".pdbqt")),
            "--cpu",
            "1",
            "--exhaustiveness",
            ("4" if arm == "pro_double" else "2"),
            "--max_evals",
            "25000",
            "--num_modes",
            "9",
            "--energy_range",
            "100",
            "--seed",
            str(seed),
            "--out",
            str(pool),
        ]
        for axis, v in zip("xyz", center):
            cmd += ["--center_" + axis, str(v), "--size_" + axis, "26"]
        if arm != "upstream":
            cmd += ["--pose_selection", "late", "--search_counters"]
        started = time.monotonic()
        with pool.with_suffix(".log").open("w") as f:
            subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT)
        wall = time.monotonic() - started
        score = float(
            next(
                l
                for l in pool.read_text().splitlines()
                if l.startswith("REMARK VINA RESULT:")
            ).split()[3]
        )
        base = dict(
            target="S6K1",
            state="3A60",
            ligand=r["cid"],
            seed=seed,
            kd_nM=r["kd_nM"],
            arm=arm,
            score=score,
            wall_s=wall,
        )
        correction = (
            run_one((args.integrator, points, pool, float(r["kd_nM"])))
            if arm == "pro"
            else []
        )
        return base, correction

    with ThreadPoolExecutor(max_workers=args.workers) as workers:
        outcomes = list(workers.map(cell, runs))
    write(args.work / "docks.csv", [d for d, _ in outcomes])
    write(args.work / "cells.csv", [x for _, corr in outcomes for x in corr])


def sensitivity(args):
    points = args.work / "independent_nodes.txt"
    nodes(points, 8192, 66173)
    rows = [
        r
        for r in csv.DictReader((args.work / "panel.csv").open())
        if r["decision"] == "included"
    ]
    inputs = [
        (
            args.integrator,
            points,
            args.work / "pools" / f"3A60_{r['cid']}_{seed}_pro.pdbqt",
            float(r["kd_nM"]),
            32,
        )
        for r in rows
        for seed in (271828, 314159)
    ]
    with ThreadPoolExecutor(max_workers=args.workers) as workers:
        data = [v for rows in workers.map(run_one, inputs) for v in rows]
    write(args.work / "sensitivity.csv", data)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=["prepare", "dock", "sensitivity"])
    p.add_argument("repo", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("--integrator", type=Path)
    p.add_argument("--pro", type=Path)
    p.add_argument("--native", type=Path)
    p.add_argument("--workers", type=int, default=1)
    args = p.parse_args()
    {"prepare": prepare, "dock": dock, "sensitivity": sensitivity}[args.action](args)
