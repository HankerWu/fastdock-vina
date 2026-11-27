"""Freeze the source graph, two preparation geometries, paired FXa searches, and within-assay rank."""

import argparse
import csv
import hashlib
from itertools import combinations, product
import re
from pathlib import Path
from statistics import median, mean
import subprocess
import sys
import tarfile
from time import perf_counter

from scipy.stats import spearmanr, kendalltau

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
sys.path.insert(0, str(BENCH.parent / "tools"))
from score_physical_pool import mapped_ligand_bonds

KEY = "CHEMBL369480"
REC = ("1FJS", "2P16")
SEED = (271828, 314159)


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def save(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def signature(path):
    return [
        l
        for l in path.read_text().splitlines()
        if l.startswith(("REMARK SMILES", "REMARK H PARENT", "BRANCH"))
    ]


def parsed(binary, ligand):
    edges = {"bond": set(), "pair": set()}
    for line in subprocess.check_output([binary, ligand], text=True).splitlines():
        tag, a, b = line.split()
        edges[tag].add((int(a), int(b)))
    return edges


def score(binary, receptor, ligand, pose, graph):
    cmd = [
        binary,
        "--receptor",
        receptor,
        "--ligand",
        ligand,
        "--pose_coordinates",
        pose,
        "--score_only",
        "--cpu",
        "1",
        "--seed",
        "214039",
    ]
    if graph:
        cmd += ["--ligand_bonds", graph]
    log = subprocess.check_output(cmd, text=True).split(
        "Estimated Free Energy of Binding"
    )[-1]
    vals = list(map(float, re.findall(r":\s+(-?\d+\.\d+)\s+\(kcal/mol\)", log)))
    if len(vals) != 10:
        raise ValueError(log)
    return vals


def label(edges):
    return " ".join(f"{a}-{b}" for a, b in sorted(edges))


def main(a):
    work = a.work
    work.mkdir(exist_ok=True, parents=True)
    for archive, members in [
        (BENCH / "frozen_core/affinity/fxa_frozen.tar.xz", None),
        (
            BENCH / "amide_search/prepared-fixtures.tar.gz",
            ["fdv-fxa/box.txt", "fdv-fxa1/box.txt"],
        ),
        (
            BENCH / "frozen_core/strain_transfer_audit/redock_raw.tar.xz",
            ["original.pdbqt", "random_coords.pdbqt"],
        ),
    ]:
        with tarfile.open(archive) as tar:
            tar.extractall(
                work,
                members=[tar.getmember(x) for x in members] if members else None,
                filter="data",
            )
    original = work / "original.pdbqt"
    clean = work / "random_coords.pdbqt"
    if signature(original) != signature(clean):
        raise ValueError("input identity mapping changed")
    bonds = mapped_ligand_bonds(clean.read_text())
    graph = work / "complete.bonds"
    graph.write_text("".join(f"{x} {y}\n" for x, y in bonds))
    audits = []
    for row in read(BENCH / "affinity_partition/fxa_series.csv"):
        if row["included"] != "1":
            continue
        key = row["ligand"]
        ligand = work / "1FJS" / (key + ".pdbqt")
        old = parsed(a.old_graph, ligand)
        new = parsed(a.fixed_graph, ligand)
        if key != KEY and old != new:
            raise ValueError("Another ligand changed; rerun it: " + key)
        if key != KEY:
            ref = mapped_ligand_bonds(ligand.read_text())
        else:
            ref = bonds
        heavy = {
            int(l[6:11])
            for l in ligand.read_text().splitlines()
            if l.startswith(("ATOM  ", "HETATM")) and l.split()[-1] not in ("H", "HD")
        }
        truth = {(i, j) for i, j in ref if i in heavy and j in heavy}
        for arm, edges in (("old", old), ("corrected", new)):
            audits.append(
                dict(
                    ligand=key,
                    arm=arm,
                    heavy_bonds=len(edges["bond"]),
                    missing=label(truth - edges["bond"]),
                    extra=label(edges["bond"] - truth),
                    pairs=len(edges["pair"]),
                    removed_pairs=label(old["pair"] - edges["pair"]),
                )
            )
    save(work / "graphs.csv", audits)
    allcells = read(BENCH / "terminal_rotor_score/results/poses.csv")
    base = [r for r in allcells if r["target"] == "FXa"]
    cells = []
    terms = []
    for geom, rec, seed in product(("original", "random_coords"), REC, SEED):
        folder = work / rec
        ligand = work / (geom + ".pdbqt")
        receptor = folder / (rec + ".pdbqt")
        cfg = work / ("fdv-fxa" if rec == "2P16" else "fdv-fxa1") / "box.txt"
        ref = folder / f"{KEY}_{seed}.pdbqt"
        pose = work / "frozen.pdbqt"
        pose.write_text(ref.read_text().split("ENDMDL")[0] + "ENDMDL\n")
        # All arms here use *original* fixed Cartesian poses, even for the new preparation.
        for arm, binary in (
            ("old", a.old),
            ("corrected", a.fixed),
            ("complete", a.fixed),
        ):
            key = f"{geom}_{rec}_{seed}_{arm}"
            output = work / (key + ".pdbqt")
            counter = [
                "--search_counters",
                "--pose_selection",
                "late",
                "--max_evals",
                "50000",
                "--num_modes",
                "9",
                "--energy_range",
                "100",
                "--verbosity",
                "0",
                "--exhaustiveness",
                "4",
                "--cpu",
                "1",
                "--seed",
                str(seed),
            ]
            cmd = [
                binary,
                "--receptor",
                receptor,
                "--ligand",
                ligand,
                "--config",
                cfg,
                *counter,
                "--out",
                output,
            ]
            if arm == "complete":
                cmd += ["--ligand_bonds", graph]
            tic = perf_counter()
            log = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
            elapsed = perf_counter() - tic
            (work / (key + ".log")).write_text(log)
            lines = output.read_text()
            top = float(re.search(r"REMARK VINA RESULT:\s+(\S+)", lines)[1])
            cnt = {
                k: sum(map(int, re.findall(r"\b" + k + r"=(\d+)", log)))
                for k in ("hunt_evals", "second_evals", "explicit_evals")
            }
            cells.append(
                dict(
                    geometry=geom,
                    receptor=rec,
                    seed=seed,
                    arm=arm,
                    score=top,
                    seconds=elapsed,
                    modes=lines.count("ENDMDL"),
                    sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                    **cnt,
                )
            )
            # Native fixed-pose API requires identical atom-identity and torsion-tree records.
            v = score(
                binary, receptor, ligand, pose, graph if arm == "complete" else None
            )
            terms.append(
                dict(
                    geometry=geom,
                    receptor=rec,
                    seed=seed,
                    arm=arm,
                    self=v[0],
                    inter=v[1],
                    intra=v[4],
                    raw_zero=v[1] + v[4],
                )
            )
    save(work / "cells.csv", cells)
    save(work / "fixed.csv", terms)
    stats = []
    comparisons = []
    for geom, arm, rec in product(
        ("original", "random_coords"), ("old", "corrected", "complete"), (*REC, "both")
    ):
        rs = REC if rec == "both" else (rec,)
        values = {}
        measured = {r["ligand"]: float(r["measurement"]) for r in base}
        for compound in measured:
            grid = []
            for r in rs:
                if compound == KEY:
                    vals = [
                        x["score"]
                        for x in cells
                        if (x["geometry"], x["arm"], x["receptor"]) == (geom, arm, r)
                    ]
                else:
                    vals = [
                        float(x["fixed"])
                        for x in base
                        if x["ligand"] == compound and x["receptor"] == r
                    ]
                assert len(vals) == 2
                grid.append(median(vals))
            values[compound] = mean(grid)
        pairs = [
            (i, j)
            for i, j in combinations(sorted(measured, key=measured.get), 2)
            if measured[j] >= 3 * measured[i] - 1e-12
        ]

        def credit(pair):
            d = values[pair[1]] - values[pair[0]]
            return 1 if d > 0 else 0.5 if d == 0 else 0

        incident = [p for p in pairs if KEY in p]
        for i, j in incident:
            comparisons.append(
                dict(
                    geometry=geom,
                    arm=arm,
                    receptor=rec,
                    strong=i,
                    weak=j,
                    gap=values[j] - values[i],
                    credit=credit((i, j)),
                )
            )
        stats.append(
            dict(
                geometry=geom,
                arm=arm,
                receptor=rec,
                pairs=len(pairs),
                correct=sum(map(credit, pairs)),
                incident_n=len(incident),
                incident_correct=sum(map(credit, incident)),
                spearman=spearmanr(
                    list(measured.values()), [values[k] for k in measured]
                ).statistic,
                kendall=kendalltau(
                    list(measured.values()), [values[k] for k in measured]
                ).statistic,
            )
        )
    save(work / "assay.csv", stats)
    save(work / "incident_pairs.csv", comparisons)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("work", "old", "fixed", "old_graph", "fixed_graph"):
        p.add_argument(key, type=lambda s: Path(s).resolve())
    main(p.parse_args())
