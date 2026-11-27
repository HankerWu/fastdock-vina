"""Matched native Late docking after a single same-seed change in coordinate initialization."""

import argparse
import hashlib
import re
import subprocess
import tarfile
from itertools import product
from pathlib import Path
from statistics import median
from tempfile import TemporaryDirectory
from time import perf_counter

import replay
from decompose import read, read_first
from evaluate import physical
from geometry import inspect, prepare
from summarize import write

LIGAND = "CHEMBL369480"
ARMS = ("original", "random_coords")


def recombine(out, cells):
    stored = read(replay.HERE.parent / "strain_components/scores.csv")
    eligible = read(out / "pairs.csv.xz")
    summaries = []
    for arm in ARMS:
        scores = {
            (r["compound"], r["receptor"], str(r["seed"])): r
            for r in stored
            if r["assay"] == "CHEMBL661804"
        }
        scores.update(
            {
                (LIGAND, r["receptor"], str(r["seed"])): r
                for r in cells
                if r["arm"] == arm
            }
        )
        for rec in ("1FJS", "2P16", "both"):
            pairs = [
                r
                for r in eligible
                if r["assay"] == "CHEMBL661804" and r["receptor"] == rec
            ]
            block = list(
                product(
                    ("1FJS", "2P16") if rec == "both" else [rec], ("271828", "314159")
                )
            )
            for subset in ("all", "incident", "remaining"):
                rows = [
                    r
                    for r in pairs
                    if subset == "all"
                    or (LIGAND in (r["stronger"], r["weaker"]))
                    == (subset == "incident")
                ]
                totals = dict.fromkeys(("native", "full"), 0)
                for r in rows:
                    for field in totals:
                        gap = median(
                            float(scores[r["weaker"], *k][field])
                            - float(scores[r["stronger"], *k][field])
                            for k in block
                        )
                        credit = 0.5 if abs(gap) < 1e-10 else int(gap > 0)
                        if arm == "original" and credit != float(r[field]):
                            raise ValueError(
                                "Original control changed an archived pair"
                            )
                        totals[field] += credit
                summaries.append(
                    dict(
                        arm=arm, receptor=rec, subset=subset, pairs=len(rows), **totals
                    )
                )
    write(out / "redock_summary.csv", summaries)


def main(binary, out):
    series = [
        r
        for r in read(replay.previous.BENCH / "affinity_partition/fxa_series.csv")
        if r["included"] == "1"
    ]
    meta = next(
        r
        for r in read(replay.HERE.parent / "strain_components/scores.csv")
        if r["compound"] == LIGAND
    )
    free = float(meta["free"])
    all_preps, cells = [], []
    with TemporaryDirectory() as tmp:
        work = Path(tmp)
        raw = work / "raw"
        raw.mkdir()
        for file in (
            "amide_search/prepared-fixtures.tar.gz",
            "frozen_core/affinity/fxa_frozen.tar.xz",
        ):
            with tarfile.open(replay.previous.BENCH / file) as archive:
                archive.extractall(work, filter="data")
        for r in series:
            c = r["ligand"]
            for arm in ARMS:
                start = perf_counter()
                _, quality, text = prepare(r["docked_smiles"], arm == "random_coords")
                sec = perf_counter() - start
                if (
                    arm == "original"
                    and text != (work / "2P16" / (c + ".pdbqt")).read_text()
                ):
                    raise ValueError(f"Original input mismatch: {c}")
                all_preps.append(
                    dict(compound=c, arm=arm, geometry_and_prep_s=sec, **quality)
                )
                if c == LIGAND:
                    (raw / (arm + ".pdbqt")).write_text(text)
        for receptor, seed in product(("1FJS", "2P16"), (271828, 314159)):
            folder = work / ("fdv-fxa" if receptor == "2P16" else "fdv-fxa1")
            receptor_file = work / receptor / (receptor + ".pdbqt")
            for arm in ARMS if seed == 271828 else reversed(ARMS):
                name = f"{arm}_{receptor}_{seed}"
                path = raw / (name + ".pdbqt")
                cmd = [
                    binary,
                    "--receptor",
                    receptor_file,
                    "--ligand",
                    raw / (arm + ".pdbqt"),
                    "--config",
                    folder / "box.txt",
                    "--seed",
                    str(seed),
                    "--cpu",
                    "1",
                    "--max_evals",
                    "50000",
                    "--exhaustiveness",
                    "4",
                    "--pose_selection",
                    "late",
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--verbosity",
                    "0",
                    "--search_counters",
                    "--out",
                    path,
                ]
                start = perf_counter()
                proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
                sec = perf_counter() - start
                log = proc.stdout + proc.stderr
                (raw / (name + ".log")).write_text(log)
                if (
                    arm == "original"
                    and path.read_bytes()
                    != (work / receptor / f"{LIGAND}_{seed}.pdbqt").read_bytes()
                ):
                    raise ValueError(f"Original exported pool differs: {name}")
                pose = read_first(path)
                # Distorted input geometry changes the historical engine's inferred covalent graph.
                scale = physical.evaluate_coordinates(
                    binary, raw / (arm + ".pdbqt"), receptor_file, None, 60
                )["scale"]
                inter = next(float(line.split()[-1]) for line in path.read_text().splitlines()
                             if line.startswith("REMARK INTER:"))
                # Each saved PDBQT term is rounded independently to three decimals.
                if abs(scale * (inter + pose.native_intra - pose.native_unbound) - pose.native_score) > 0.002:
                    raise ValueError(f"Normalization disagrees with saved native score: {name}")
                mol, _ = physical.decode_ligand(pose)
                bound, _, _ = replay.fixed_geometry(pose)
                costs = {}
                for k in (
                    "hunt_evals",
                    "second_evals",
                    "explicit_evals",
                    "hunt_calls",
                    "second_calls",
                    "explicit_calls",
                ):
                    hits = re.findall(r"\b" + k + r"=(\d+)", log)
                    if len(hits) != (1 if k.startswith("explicit") else 4):
                        raise ValueError(f"Incomplete accounting: {name} {k}")
                    costs[k] = sum(map(int, hits))
                cells.append(
                    dict(
                        arm=arm,
                        receptor=receptor,
                        seed=seed,
                        native=pose.native_score,
                        scale=scale,
                        full=physical.mmff_affinity_score(pose, scale, bound, free),
                        bound=bound,
                        free=free,
                        max_bond_deviation_A=inspect(mol)["max_bond_deviation_A"],
                        wall_s=sec,
                        modes=len(physical.read_poses(path)),
                        **costs,
                    )
                )
        with tarfile.open(out / "redock_raw.tar.xz", "w:xz") as archive:
            for file in sorted(raw.iterdir()):
                info = archive.gettarinfo(file, arcname=file.name)
                info.mtime = info.uid = info.gid = 0
                info.uname = info.gname = ""
                with file.open("rb") as stream:
                    archive.addfile(info, stream)
    write(out / "input_geometry.csv", all_preps)
    write(out / "redock_cells.csv", cells)
    recombine(out, cells)
    print(f"binary sha256 {hashlib.sha256(binary.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    main(args.binary.resolve(), args.out.resolve())
