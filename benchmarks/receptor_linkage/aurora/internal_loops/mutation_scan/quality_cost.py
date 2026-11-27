"""Check finite chemistry and independently time uninstrumented paired commands."""

import argparse
import csv
import os
import re
import resource
import subprocess
import time
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree
import run


def pdbqt_atoms(lines):
    return {
        int(l[6:11]): np.array([float(l[j : j + 8]) for j in (30, 38, 46)])
        for l in lines
        if l.startswith(("ATOM  ", "HETATM"))
    }


def summary(aurora, brd9, study):
    with (study / "cells.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    inputs = run.config(aurora, brd9, study)
    qa = []
    for drug, cell in inputs.items():
        types = {
            int(l[6:11]): l.split()[-1]
            for l in cell["ligand"].read_text().splitlines()
            if l.startswith(("ATOM  ", "HETATM"))
        }
        heavy = set(cell["serials"])
        baseline = pdbqt_atoms(cell["ligand"].read_text().splitlines())
        edges = [
            tuple(map(int, l.split())) for l in cell["bonds"].read_text().splitlines()
        ]
        edges = [(a, b) for a, b in edges if a in heavy and b in heavy]
        if drug == "82I":
            stereo = []
        else:
            mol = Chem.MolFromMolFile(
                str(cell["ligand"].with_suffix(".sdf")), removeHs=False
            )
            _, ids, indices, _, _ = run.get_tree(cell["ligand"])
            serial = dict(zip(indices, ids, strict=True))
            stereo = []
            for idx, tag in Chem.FindMolChiralCenters(mol, includeUnassigned=True):
                n = [
                    serial[a.GetIdx()]
                    for a in mol.GetAtomWithIdx(idx).GetNeighbors()
                    if a.GetAtomicNum() > 1
                ]
                if tag == "?" or len(n) != 3:
                    raise ValueError(f"Unresolved stereochemistry {drug} {idx}")
                stereo.append((serial[idx], n))

        def chiral(x):
            return [
                np.sign(np.linalg.det([x[n] - x[i] for n in nn])) for i, nn in stereo
            ]

        def bond_lengths(x):
            return np.array([np.linalg.norm(x[a] - x[b]) for a, b in edges])

        receptor = pdbqt_atoms(
            l
            for l in cell["receptor"].read_text().splitlines()
            if l.split() and l.split()[-1] not in ("HD", "H")
        )
        tree = cKDTree(list(receptor.values()))
        for row in (r for r in rows if r["drug"] == drug):
            pool = study / "raw" / f"{drug}_{row['seed']}_{row['arm']}.pdbqt"
            for rank, block in enumerate(pool.read_text().split("ENDMDL"), 1):
                xyz = pdbqt_atoms(block.splitlines())
                if not xyz:
                    continue
                finite = bool(np.isfinite(list(xyz.values())).all())
                if set(xyz) != set(types) or not finite:
                    raise ValueError(f"Malformed output {pool} {rank}")
                nearest = tree.query([xyz[i] for i in sorted(heavy)])[0].min()
                qa.append(
                    dict(
                        drug=drug,
                        seed=row["seed"],
                        arm=row["arm"],
                        rank=rank,
                        all_atoms_finite=int(finite),
                        max_heavy_bond_drift_A=abs(
                            bond_lengths(xyz) - bond_lengths(baseline)
                        ).max(),
                        inverted_centres=sum(
                            a != b
                            for a, b in zip(chiral(xyz), chiral(baseline), strict=True)
                        ),
                        centres_checked=len(stereo),
                        closest_receptor_heavy_A=nearest,
                    )
                )
    run.replay.geom.write_csv(study, "qa.csv", qa)
    return inputs, rows


def timing(inputs, rows, native, scan, study):
    timings = []
    raw = study / "cost"
    raw.mkdir(exist_ok=False)
    env = {
        k: v for k, v in os.environ.items() if k not in ("FD_SCAN_TRACE", "FD_SCAN_IID")
    }
    for row in (r for r in rows if r["arm"] == "native"):
        drug, seed = row["drug"], int(row["seed"])
        for repeat in (0, 1):
            arms = ("native", "scan") if repeat == 0 else ("scan", "native")
            for arm in arms:
                output = raw / f"{drug}_{seed}_{arm}_{repeat}.pdbqt"
                expected = study / "raw" / f"{drug}_{seed}_{arm}.pdbqt"
                before = resource.getrusage(resource.RUSAGE_CHILDREN)
                start = time.perf_counter()
                with output.with_suffix(".log").open("w") as f:
                    subprocess.run(
                        run.command(
                            native if arm == "native" else scan,
                            inputs[drug],
                            seed,
                            output,
                        ),
                        stdout=f,
                        stderr=subprocess.STDOUT,
                        env=env,
                        check=True,
                    )
                wall = time.perf_counter() - start
                after = resource.getrusage(resource.RUSAGE_CHILDREN)
                if output.read_bytes() != expected.read_bytes():
                    raise ValueError(f"Timing replay diverged at {output.name}")
                log = output.with_suffix(".log").read_text()
                evals = sum(
                    map(int, re.findall(r"\b(?:hunt|second|explicit)_evals=(\d+)", log))
                )
                timings.append(
                    dict(
                        drug=drug,
                        seed=seed,
                        arm=arm,
                        repeat=repeat,
                        evals=evals,
                        wall_s=wall,
                        cpu_s=after.ru_utime
                        + after.ru_stime
                        - before.ru_utime
                        - before.ru_stime,
                        pool_sha256=run.replay.digest(output),
                    )
                )
    run.replay.geom.write_csv(study, "timing.csv", timings)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("aurora", "brd9", "study", "native", "scan"):
        p.add_argument(arg, type=Path)
    args = p.parse_args()
    cells, rows = summary(args.aurora, args.brd9, args.study)
    timing(cells, rows, args.native.resolve(), args.scan.resolve(), args.study)
