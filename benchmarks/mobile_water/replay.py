"""Replay archived dry candidate pools. See README for inputs and environment."""

import argparse
import csv
import itertools
import json
import resource
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import kendalltau, spearmanr

from mobile_water import Site

METHODS = ("dry", "point", "mobile", "order17")
LANDMARKS = ("57", "60A", "60D", "99", "174", "189", "195", "215", "216", "227", "228")


def records(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def save(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def xyz(line):
    return [float(line[i : i + 8]) for i in (30, 38, 46)]


def pdb_lines(path):
    return [
        l
        for l in path.read_text().splitlines()
        if l.startswith(("ATOM  ", "HETATM"))
        and l[16] in " A"
        and l[76:78].strip() not in ("H", "D")
    ]


def thrombin_site(source, receiver):
    def frame(rows):
        c = {
            l[22:27].strip(): xyz(l)
            for l in rows
            if l.startswith("ATOM") and l[21] == "H" and l[12:16] == " CA "
        }
        return np.array([c[key] for key in LANDMARKS])

    a, b = frame(source), frame(receiver)
    u, _, v = np.linalg.svd((a - a.mean(0)).T @ (b - b.mean(0)))
    u[:, -1] *= np.linalg.det(u @ v)
    (w,) = [l for l in source if l[17:20] == "HOH" and l[21:26] == "H2006"]
    return [("3D49:H2006", (xyz(w) - a.mean(0)) @ u @ v + b.mean(0), float(w[60:66]))]


def cross_sites(rows, protein, center, chain="A"):
    polar = np.array(
        [
            xyz(l)
            for l in rows
            if l.startswith("ATOM")
            and l[21] == chain
            and l[76:78].strip() in ("N", "O")
        ]
    )
    metals = np.array(
        [
            xyz(l)
            for l in rows
            if l[21] == chain
            and l[76:78].strip() in ("ZN", "MG", "MN", "FE", "CA", "CO", "NI", "CU")
        ]
    )
    selected = []
    for l in rows:
        if l[17:20] != "HOH" or l[21] != chain or float(l[54:60]) < 0.5:
            continue
        site = np.array(xyz(l))
        b_iso = float(l[60:66])
        if b_iso <= 0 or max(abs(site - center)) > 12 or protein.query(site)[0] < 2.3:
            continue
        distances = np.linalg.norm(polar - site, axis=1)
        if np.count_nonzero((distances >= 2.4) & (distances <= 3.5)) < 2:
            continue
        if metals.size and np.linalg.norm(metals - site, axis=1).min() < 3.0:
            continue
        selected.append((l[21:27].strip(), site, b_iso))
    return selected


def models(path):
    block = path.read_text().split("ENDMDL")
    if block[-1].strip():
        raise ValueError(f"Incomplete pool {path}")
    for body in block[:-1]:
        lines = body.strip().splitlines()
        model = int(lines[0].split()[1])
        (score,) = [
            float(l.split()[3]) for l in lines if l.startswith("REMARK VINA RESULT:")
        ]
        position = []
        for line in lines:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            ad_type = line.split()[-1]
            if ad_type not in ("H", "HD", "HS", "W") and not (
                ad_type.startswith("G") and ad_type[1:].isdigit()
            ):
                position.append(xyz(line))
        if not position:
            raise ValueError(f"Missing heavy atoms in {path}, model {model}")
        yield model, score, cKDTree(position)


def input_data(args):
    data = {}
    sources = {}
    cross = args.repo / "benchmarks" / "paired_crossdocking" / "results"
    throm = args.repo / "benchmarks" / "thrombin_water" / "expected"
    for r in records(throm / "poses.csv"):
        rec, lig, seed = r["rec"], r["ligand"], int(r["seed"])
        key = ("thrombin", 50000, rec, lig, seed)
        path = args.throm / "results" / f"{rec}_{lig}_{seed}.pdbqt"
        data.setdefault(key, (path, {}))[1][int(r["model"])] = (
            r["rmsd"],
            float(r["native"]),
            float(r["exclusion"]),
        )
    for budget in (3000, 12000):
        for r in records(cross / f"poses_{budget}.csv"):
            if r["arm"] != "native" or r["kind"] != "cross":
                continue
            rec, lig, seed = r["receiver"], r["ligand"], int(r["seed"])
            key = (r["target"], budget, rec, lig, seed)
            path = args.cross / f"runs_{budget}" / f"{lig}_{rec}_{seed}_native.pdbqt"
            data.setdefault(key, (path, {}))[1][int(r["rank"])] = (
                r["rmsd"],
                float(r["score"]),
                None,
            )
    if {k[0] for k in data} != {"thrombin", "PDE4D", "ESR1", "PPARG"} or len(
        data
    ) != 102:
        raise ValueError("Incomplete registered panel")
    metadata = json.loads((args.cross / "metadata.json").read_text())
    for target, _, rec, _, _ in data:
        if rec in sources:
            continue
        dev = target == "thrombin"
        base = args.throm if dev else args.cross
        rows = pdb_lines(base / "raw" / f"{rec}.pdb")
        prot = cKDTree(
            [
                xyz(l)
                for l in rows
                if l.startswith("ATOM") and l[21] in ("HL" if dev else "A")
            ]
        )
        if dev:
            sites = thrombin_site(pdb_lines(base / "raw" / "3D49.pdb"), rows)
        else:
            sites = cross_sites(rows, prot, np.array(metadata[rec]["center"]))
        sources[rec] = (target, prot, sites)
    return data, sources, records(throm / "series.csv")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("repo", "throm", "cross", "out"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    cpu = time.process_time()
    data, sources, assay = input_data(args)
    site_rows, prepared, costs = [], {}, defaultdict(float)
    for rec, (target, protein, sites) in sources.items():
        for label, coords, b in sites:
            site_rows.append(
                {
                    "target": target,
                    "receiver": rec,
                    "site": label,
                    "x": coords[0],
                    "y": coords[1],
                    "z": coords[2],
                    "b_iso": b,
                }
            )
        prepared[rec] = {}
        for order in (11, 17):
            then = time.perf_counter()
            prepared[rec][order] = [
                Site.prepare(c, b, protein, order) for _, c, b in sites
            ]
            costs[(target, f"receiver_order{order}")] += time.perf_counter() - then
    save(args.out / "sites.csv", site_rows)
    poses, cells, scores, timing = [], [], {}, []
    for key, (path, labels) in sorted(data.items()):
        target, budget, receiver, ligand, seed = key
        meta = {
            "target": target,
            "budget": budget,
            "receiver": receiver,
            "ligand": ligand,
            "seed": seed,
        }
        pool = []
        for model, dry, tree in models(path):
            rmsd, archived, published_point = labels[model]
            if abs(dry - archived) > 1e-8:
                raise ValueError(f"Native score mismatch in {path}, model {model}")
            then = time.perf_counter()
            primary = [s.evaluate(tree) for s in prepared[receiver][11]]
            costs[(target, "evaluate_order11")] += time.perf_counter() - then
            then = time.perf_counter()
            higher = [s.evaluate(tree) for s in prepared[receiver][17]]
            costs[(target, "evaluate_order17")] += time.perf_counter() - then
            point = dry + sum(v[0] for v in primary)
            if published_point is not None and abs(point - published_point) > 1e-7:
                raise ValueError(
                    f"Historical rigid-site ablation mismatch in {path}: {point} {published_point}"
                )
            pool.append(
                dict(
                    **meta,
                    model=model,
                    rmsd=float(rmsd) if rmsd else None,
                    dry=dry,
                    point=point,
                    mobile=dry + sum(v[1] for v in primary),
                    order17=dry + sum(v[1] for v in higher),
                )
            )
        if {p["model"] for p in pool} != set(labels):
            raise ValueError(f"Incomplete candidate correspondence in {path}")
        poses.extend(pool)
        for arm in METHODS:
            ranked = sorted(pool, key=lambda p: (p[arm], p["model"]))
            first = ranked[0]
            scores[(*key, arm)] = first[arm]
            if first["rmsd"] is not None:
                cells.append(
                    dict(
                        **meta,
                        arm=arm,
                        model=first["model"],
                        rmsd=first["rmsd"],
                        top1=int(first["rmsd"] <= 2),
                        top5=int(any(p["rmsd"] <= 2 for p in ranked[:5])),
                        oracle=int(any(p["rmsd"] <= 2 for p in ranked)),
                    )
                )
    save(args.out / "poses.csv", poses)
    save(args.out / "cells.csv", cells)
    for (t, stage), seconds in sorted(costs.items()):
        timing.append({"target": t, "stage": stage, "seconds": seconds})
    save(args.out / "timing.csv", timing)
    correlations, pairs = [], []
    ass = {r["id"]: float(r["Kd_uM_derived_at_298.15K"]) for r in assay}
    gaps = {}
    for arm in METHODS:
        for rec, seed in sorted({(k[2], k[4]) for k in data if k[0] == "thrombin"}):
            s = {i: scores[("thrombin", 50000, rec, i, seed, arm)] for i in ass}
            ordered = sorted(ass, key=ass.get)
            correlations.append(
                {
                    "arm": arm,
                    "receiver": rec,
                    "seed": seed,
                    "spearman": spearmanr(
                        [ass[i] for i in ordered], [s[i] for i in ordered]
                    ).statistic,
                    "kendall": kendalltau(
                        [ass[i] for i in ordered], [s[i] for i in ordered]
                    ).statistic,
                }
            )
            for a, b in itertools.combinations(ordered, 2):
                if ass[b] / ass[a] >= 3:
                    gaps.setdefault((arm, a, b), []).append(s[b] - s[a])
    for (arm, a, b), g in gaps.items():
        med = float(np.median(g))
        pairs.append(
            {
                "arm": arm,
                "stronger": a,
                "weaker": b,
                "median_gap": med,
                "correct": 1.0 if med > 0 else (0.5 if med == 0 else 0.0),
            }
        )
    save(args.out / "assay_pairs.csv", pairs)
    save(args.out / "assay_correlations.csv", correlations)
    elapsed = {
        "cpu_s": time.process_time() - cpu,
        "wall_s": time.perf_counter() - start,
        "peak_resident_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (args.out / "runtime.json").write_text(json.dumps(elapsed, indent=2) + "\n")
    print(
        "poses",
        len(poses),
        "cells",
        len(cells),
        "peak_resident_kib",
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    )
    for target in sorted({row["target"] for row in poses}):
        for arm in METHODS:
            group = [r for r in cells if r["target"] == target and r["arm"] == arm]
            print(
                target,
                arm,
                len(group),
                *(sum(r[k] for r in group) for k in ("top1", "top5", "oracle")),
            )
    print(
        "ITC eligible pair credits",
        [(m, sum(r["correct"] for r in pairs if r["arm"] == m)) for m in METHODS],
    )


if __name__ == "__main__":
    main()
