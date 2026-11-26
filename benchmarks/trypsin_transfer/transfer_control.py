"""Test the fixed full-protein OBC score with a second amidinium target and spatial truncation."""

import argparse
import itertools
import math
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from statistics import median

import numpy as np
from data import read, write, ROOT, SEEDS, RT, RT_FXA, LIGANDS
from select_zinc_pose import read_models
from polar_endpoint import PolarEndpoint, atoms

ARMS = ("native", "shell")


def restricted(lines, radius):
    rows = [l for l in lines if l.startswith(("ATOM  ", "HETATM"))]
    xyz = np.array([[float(l[k : k + 8]) for k in (30, 38, 46)] for l in rows])
    select = [
        i
        for i, l in enumerate(rows)
        if l[21] == "A"
        and (
            (l[17:20] == "SER" and int(l[22:26]) == 195 and l[12:16].strip() == "OG")
            or (l[17:20] == "ASP" and int(l[22:26]) == 189 and l[12:16].strip() == "CG")
        )
    ]
    assert len(select) == 2
    dist = np.linalg.norm(xyz[:, None] - xyz[select], axis=-1).min(axis=1)
    residues = {line[17:27] for line, d in zip(rows, dist, strict=True) if d < radius}
    return [l for l in rows if l[17:27] in residues]


def job(task):
    target, rec, name, folder, native_dir = task
    protein_file = native_dir / f"{rec}.pdbqt"
    lines = protein_file.read_text().splitlines()
    contacts = (
        [("full", lines)]
        if target == "fxa"
        else [
            (k, restricted(lines, v))
            for k, v in (("6A", 6.0), ("10A", 10.0), ("16A", 16.0))
        ]
    )
    ligand = atoms((native_dir / (name + ".pdbqt")).read_text().splitlines())
    results = []
    for method, subset in contacts:
        rec_atoms = atoms(subset, receptor=True)
        calc = PolarEndpoint(rec_atoms, ligand)
        for arm, seed in itertools.product(ARMS, SEEDS):
            path = folder / f"{rec}_{name}_{seed}_{arm}.pdbqt"
            if target == "fxa" and arm == "native" and seed in (271828, 314159):
                path = native_dir / f"{name}_{seed}.pdbqt"
            top = read_models(path.read_text().splitlines())[1]
            native = float(
                next(l for l in top if l.startswith("REMARK VINA RESULT")).split()[3]
            )
            values = calc.evaluate(atoms(top))
            results.append(
                dict(
                    target=target,
                    receptor=rec,
                    ligand=name,
                    seed=seed,
                    arm=arm,
                    region=method,
                    receptor_atoms=len(rec_atoms.xyz),
                    vina=native,
                    **values,
                )
            )
    return results


def gaps(data, fxa):
    sizes = {r["ligand"]: (r["heavy"], r["charge"]) for r in fxa}
    labels = {
        "fxa": {r["ligand"]: float(r["ki_nm"]) for r in fxa},
        "trypsin": {
            r["ligand"]: float(r["kd_uM"])
            for r in read(ROOT / "benchmarks/directional_ion_pair/results/series.csv")
        },
    }
    keys = ("target", "arm", "region", "receptor", "seed")
    blocks = defaultdict(dict)
    for r in data:
        blocks[tuple(r[k] for k in keys)][r["ligand"]] = r
    accumulated = defaultdict(list)
    for (target, arm, region, rec, seed), group in blocks.items():
        kds = labels[target]
        for a, b in itertools.combinations(sorted(group, key=kds.get), 2):
            for method in ("vina", "polar_transfer", "vina_plus_polar"):

                def E(l):
                    return (
                        float(group[l]["vina"]) + float(group[l]["polar_transfer"])
                        if method == "vina_plus_polar"
                        else float(group[l][method])
                    )

                for subset in (rec, "both"):
                    accumulated[target, arm, region, subset, method, a, b].append(
                        round(E(b) - E(a), 9)
                    )
    edges = []
    summary = defaultdict(list)
    for (target, arm, region, rec, method, a, b), block in accumulated.items():
        ratio = labels[target][b] / labels[target][a]
        med = median(block)
        subsets = ["all"]
        if target == "fxa" and sizes[a] == sizes[b]:
            subsets.append("same_size_charge")
        if target == "trypsin" and "CONH2" not in (a, b):
            subsets.append("unambiguous")
        edges.append(
            dict(
                target=target,
                arm=arm,
                region=region,
                receptor=rec,
                method=method,
                strong=a,
                weak=b,
                ratio=ratio,
                predicted=med,
            )
        )
        for spec in subsets:
            if ratio >= 3:
                summary[target, arm, region, rec, method, spec].append(
                    (med, (RT if target == "trypsin" else RT_FXA) * math.log(ratio))
                )
    return edges, [
        dict(
            target=t,
            arm=a,
            region=r,
            receptor=c,
            method=m,
            subset=s,
            n=len(vals),
            correct=sum((x > 0) + 0.5 * (x == 0) for x, y in vals),
            mae_kcal=sum(abs(x - y) for x, y in vals) / len(vals),
        )
        for (t, a, r, c, m, s), vals in summary.items()
    ]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("ionpair", "fjs", "p16", "results"):
        p.add_argument(name, type=Path)
    p.add_argument("--jobs", type=int, default=3)
    a = p.parse_args()
    fxa = [
        r
        for r in read(ROOT / "benchmarks/affinity_partition/fxa_series.csv")
        if r["included"] == "1"
    ]
    jobs = [
        ("fxa", rec, x["ligand"], a.ionpair / "fxa", folder)
        for rec, folder in (("1FJS", a.fjs), ("2P16", a.p16))
        for x in fxa
    ]
    jobs.extend(
        ("trypsin", rec, name, a.ionpair / "trypsin", a.ionpair / "trypsin")
        for rec, name in itertools.product(("1BTY", "3PTB"), LIGANDS)
    )
    with ProcessPoolExecutor(a.jobs) as exe:
        results = list(itertools.chain.from_iterable(exe.map(job, jobs)))
    write(a.results / "regions.csv", results)
    cells = {
        tuple(r[k] for k in ("receptor", "seed", "arm", "ligand")): r
        for r in read(a.results / "cells.csv")
    }
    for r in read(a.results / "classical.csv"):
        r.update(
            target="trypsin",
            region="full",
            vina=cells[tuple(r[k] for k in ("receptor", "seed", "arm", "ligand"))][
                "vina"
            ],
        )
        results.append(r)
    pairs, summary = gaps(results, fxa)
    write(a.results / "region_summary.csv", summary)
    write(
        a.results / "trypsin_region_edges.csv",
        [r for r in pairs if r["target"] == "trypsin"],
    )
