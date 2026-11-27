"""Read archived PDBQT predictions and evaluate the one paired human Kd comparison."""

import argparse
import csv
import itertools
import json
import math
import sys
from pathlib import Path
from statistics import median

import numpy as np
from scipy.spatial.distance import cdist


def save(path, rows):
    with path.open("w") as f:
        w = csv.DictWriter(f, rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def geometry(repo, work, out, runs):
    sys.path.insert(0, str(repo / "benchmarks/paired_crossdocking"))
    from run import alignment

    sys.path.insert(0, str(repo / "benchmarks/macrocycle_labels"))
    from geometry import coord_blocks

    meta = {
        i: json.loads((work / i / "preparation.json").read_text())
        for i in ("1FKJ", "1FKF", "6VCU")
    }
    receptor = {
        i: next(coord_blocks(work / i / "receptor.pdbqt")) for i in ("1FKJ", "1FKF")
    }
    ligands = {
        i: next(coord_blocks(work / i / "ligand.pdbqt")) for i in ("1FKJ", "6VCU")
    }
    all_modes = []
    reference_rows = {}
    for r in runs:
        q, receiver = r["query"], r["receiver"]
        # The fallback receptor also contains FK506, with exactly the same prepared input/serials.
        # Its own ligand gives a direct cognate reference in that cell.
        source = receiver if q == "1FKJ" else q
        rot, offset, ca_rmsd, nca = alignment(meta[source], meta[receiver])
        ref = np.array(meta[source]["mapped"]) @ rot + offset
        if meta[source]["serial"] != meta[q]["serial"]:
            raise ValueError("Cognate inputs have inconsistent atom serials")
        ends = [j for j, (ty, _) in ligands[q].items() if ty == "CG0"]
        if len(ends) != 2:
            raise ValueError("Expected a single recorded ring cut")
        rec = np.asarray(
            [pos for ty, pos in receptor[receiver].values() if ty not in ("H", "HD")]
        )
        reference_distances = cdist(ref[0], rec)
        reference_rows[q, receiver] = dict(
            query=q,
            receiver=receiver,
            reference=source,
            reference_alternatives=len(ref),
            ca_rmsd=ca_rmsd,
            n_ca=nca,
            closest_heavy_A=reference_distances.min(),
            heavy_pairs_below_2A=int((reference_distances < 2).sum()),
        )
        file = out / "raw" / (r["key"] + ".pdbqt")
        score = [
            float(s.split()[3])
            for s in file.read_text().splitlines()
            if s.startswith("REMARK VINA RESULT:")
        ]
        for j, pose in enumerate(coord_blocks(file)):
            xyz = np.asarray([pose[v][1] for v in meta[q]["serial"]])
            d = cdist(xyz, rec)
            rmsd = float(
                np.sqrt(np.mean(np.sum((xyz - ref) ** 2, axis=-1), axis=-1)).min()
            )
            bond = np.linalg.norm(pose[ends[0]][1] - pose[ends[1]][1])
            all_modes.append(
                dict(
                    key=r["key"],
                    arm=r["arm"],
                    receiver=receiver,
                    query=q,
                    seed=r["seed"],
                    cap=r["cap"],
                    rank=j + 1,
                    score=score[j],
                    rmsd=rmsd,
                    closure_A=bond,
                    closest_heavy_A=float(d.min()),
                    heavy_pairs_below_2A=int((d < 2).sum()),
                    cognate=source == receiver,
                    ca_rmsd=ca_rmsd,
                    nca=nca,
                )
            )
    save(out / "reference_geometry.csv", list(reference_rows.values()))
    save(out / "modes.csv", all_modes)
    return all_modes


def main(repo, work, out):
    with (out / "runs.csv").open() as f:
        runs = list(csv.DictReader(f))
    modes = geometry(repo, work, out, runs)
    ordered = []
    for row in runs:
        pool = [m for m in modes if m["key"] == row["key"]]
        allowed = [m for m in pool if m["closure_A"] <= 3]
        before = pool[0]
        after = allowed[0] if allowed else before
        ordered.append(
            dict(
                key=row["key"],
                arm=row["arm"],
                pool=len(pool),
                within_limit=len(allowed),
                all_over_limit=not allowed,
                changed=before["rank"] != after["rank"],
                original_rank_after=after["rank"],
                before_rmsd=before["rmsd"],
                after_rmsd=after["rmsd"],
                before_closure_A=before["closure_A"],
                after_closure_A=after["closure_A"],
                before_success=before["rmsd"] <= 2,
                after_success=after["rmsd"] <= 2,
            )
        )
    save(out / "fixed_closure_order.csv", ordered)
    arms = sorted({r["arm"] for r in runs})
    exp = 0.00198720425864083 * 298.15 * math.log(119.4 / 2.7)
    best = {
        (r["arm"], r["receiver"], r["cap"], r["seed"], r["query"]): float(
            r["top_score"]
        )
        for r in runs
    }
    gaps = []
    for arm, rec, cap in itertools.product(arms, ("1FKJ", "1FKF"), ("3000", "12000")):
        seed_results = []
        for seed in ("740317", "191009", "median"):
            if seed == "median":
                strong, weak = (
                    median(row[k] for row in seed_results)
                    for k in ("fk506_score", "apx879_score")
                )
            else:
                strong = best[arm, rec, cap, seed, "1FKJ"]
                weak = best[arm, rec, cap, seed, "6VCU"]
            gap = weak - strong
            row = dict(
                arm=arm,
                receiver=rec,
                cap=cap,
                seed=seed,
                fk506_score=strong,
                apx879_score=weak,
                weak_minus_strong_score=round(gap, 6),
                measured_fold=119.4 / 2.7,
                measured_ddG_kcal_mol=exp,
                pair_credit=(gap > 0) + 0.5 * (gap == 0),
            )
            if seed != "median":
                seed_results.append(row)
            gaps.append(row)
    save(out / "kd.csv", gaps)
    summaries = []
    for arm in arms:
        obs = [x for x in modes if x["arm"] == arm]
        cells = [r for r in runs if r["arm"] == arm]
        for rname in ("1FKJ", "6VCU", "both"):
            group = [r for r in cells if rname == "both" or r["query"] == rname]
            pools = [[m for m in obs if m["key"] == r["key"]] for r in group]
            first = [p[0] for p in pools]
            oracle = sum(min(m["rmsd"] for m in p) <= 2 for p in pools)
            hits = sum(m["rmsd"] <= 2 for m in first)
            summaries.append(
                dict(
                    arm=arm,
                    query=rname,
                    cells=len(group),
                    modes=sum(map(len, pools)),
                    top1=hits,
                    top5=sum(min(m["rmsd"] for m in p[:5]) <= 2 for p in pools),
                    oracle=oracle,
                    top1_if_oracle=hits / oracle if oracle else "",
                    median_top_rmsd=median(m["rmsd"] for m in first),
                    median_oracle_rmsd=median(min(m["rmsd"] for m in p) for p in pools),
                    open_top1_above_3A=sum(m["closure_A"] > 3 for m in first),
                    open_all_above_3A=sum(m["closure_A"] > 3 for p in pools for m in p),
                    tops_clashing_2A=sum(m["closest_heavy_A"] < 2 for m in first),
                    joint=sum(int(r["joint"]) for r in group),
                    wall_s=sum(int(r["wall_ns"]) for r in group) / 1e9,
                )
            )
    save(out / "geometry_cost.csv", summaries)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("repo", "work", "out"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    main(a.repo.resolve(), a.work, a.out)
