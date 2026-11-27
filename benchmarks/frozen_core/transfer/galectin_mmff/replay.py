"""Apply published MMFF replacements to fixed native lectin pools; no fitted settings."""

import argparse
import csv
import hashlib
import itertools
import json
import lzma
import sys
import time
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from statistics import mean, median

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from scipy.stats import kendalltau, spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "benchmarks/frozen_core/strain"))
from pose import physical, score_pose

GAL = ROOT / "benchmarks/galectin_affinity"
TARGETS = {"Gal1": ("1GZW", "4Q26"), "Gal3": ("1KJL", "3T1L")}
SEEDS = (271828, 314159, 161803, 141421)
ARMS = ("native", "full", "prepared")


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    stream = lzma.open if path.suffix == ".xz" else open
    with stream(path, "wt", newline="") as f:
        out = csv.DictWriter(f, fieldnames=list(rows[0]))
        out.writeheader()
        out.writerows(rows)


def score(work, binary):
    geometries = {
        (r["pdb"], int(r["compound"]), int(r["seed"]), int(r["mode"])): r
        for r in read(GAL / "results/poses.csv")
        if r["arm"] == "native"
    }
    results = []
    for compound in range(5, 18):
        input_path = work / "prepared" / f"{compound}.pdbqt"
        input_pose = physical.Pose(1, 0, 1, input_path.read_text(), 0, 0)
        mol, identity = physical.decode_ligand(input_pose)
        free, nfree = physical.gas_reference(mol)
        scale = physical.evaluate_coordinates(
            binary, input_path, work / "prepared/1GZW.pdbqt", None, 90
        )["scale"]
        prepared = score_pose(input_pose, identity, free, nfree, scale)["fixed_mmff"]
        for pdb, seed in itertools.product(sum(TARGETS.values(), ()), SEEDS):
            pool = work / "pools" / f"{pdb}_{compound}_{seed}_native.pdbqt"
            digest = hashlib.sha256(pool.read_bytes()).hexdigest()
            assert (
                digest
                == json.loads(pool.with_suffix(".json").read_text())["out_sha256"]
            )
            poses = physical.read_poses(pool)
            physical.fixed_coordinate_inputs(poses, input_path, None)
            inter = [
                float(line.split()[-1])
                for line in pool.read_text().splitlines()
                if line.startswith("REMARK INTER:")
            ]
            for p, intermolecular in zip(poses, inter, strict=True):
                ref = geometries.pop((pdb, compound, seed, p.model))
                assert (
                    p.native_score == float(ref["score"]) and p.native_rank == p.model
                )
                r = score_pose(p, identity, free, nfree, scale)
                bound = r["fixed_mmff"]
                full = r["full_unit"]
                direct = scale * (intermolecular + bound - free)
                if abs(direct - full) > 0.002:
                    raise ValueError(
                        "Independent rounded intermolecular reconstruction diverges"
                    )
                results.append(
                    dict(
                        compound=compound,
                        pdb=pdb,
                        seed=seed,
                        model=p.model,
                        core_rmsd=float(ref["core_rmsd"]),
                        scale=scale,
                        nfree=nfree,
                        bound=bound,
                        free=free,
                        input_energy=prepared,
                        intra=p.native_intra,
                        unbound=p.native_unbound,
                        native=p.native_score,
                        full=full,
                        prepared=physical.mmff_affinity_score(
                            p, scale, bound, prepared
                        ),
                        inter_rounding_error=direct - full,
                    )
                )
    if geometries:
        raise ValueError("Archived native pool has unscored modes")
    return results


def evaluate(rows, work, out, arms=ARMS, write_table=write):
    label = {int(r["compound"]): r for r in read(work / "panel.csv")}
    ids = sorted(label)
    formula = {
        i: rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(label[i]["smiles"]))
        for i in ids
    }
    kd_table = {int(r["compound"]): r for r in read(GAL / "table1.csv")}
    pools = defaultdict(list)
    for row in rows:
        pools[row["pdb"], row["compound"], row["seed"]].append(row)
    assert len(pools) == 208
    geometry, pairs, stats = [], [], []
    for target, states in TARGETS.items():
        kd = {i: Decimal(kd_table[i][target + "_Kd_uM"]) for i in ids}
        assert all(k == Decimal(label[i][target + "_Kd_uM"]) for i, k in kd.items())
        for arm in arms:
            per_state = {}
            for pdb in states:
                top = top5 = oracle = 0
                for i, seed in itertools.product(ids, SEEDS):
                    ranked = sorted(
                        pools[pdb, i, seed], key=lambda r: (r[arm], r["model"])
                    )
                    top += ranked[0]["core_rmsd"] <= 2
                    top5 += any(r["core_rmsd"] <= 2 for r in ranked[:5])
                    oracle += any(r["core_rmsd"] <= 2 for r in ranked)
                geometry.append(
                    dict(
                        target=target,
                        pdb=pdb,
                        arm=arm,
                        pools=len(ids) * len(SEEDS),
                        top1=top,
                        top5=top5,
                        oracle=oracle,
                    )
                )
                per_state[pdb] = {
                    i: median(pools[pdb, i, s][0][arm] for s in SEEDS) for i in ids
                }
            per_state["mean"] = {i: mean(per_state[p][i] for p in states) for i in ids}
            for pdb, scores in per_state.items():
                edges = []
                for i, j in itertools.permutations(ids, 2):
                    if kd[j] < 3 * kd[i]:
                        continue
                    diff = scores[j] - scores[i]
                    edge = dict(
                        target=target,
                        pdb=pdb,
                        arm=arm,
                        stronger=i,
                        weaker=j,
                        fold=float(kd[j] / kd[i]),
                        gap=diff,
                        triazoles=17 not in (i, j),
                        same_formula=formula[i] == formula[j],
                        credit=0.5 if abs(diff) < 1e-9 else int(diff > 0),
                    )
                    edges.append(edge)
                    pairs.append(edge)
                for subset in ("all", "triazoles", "same_formula"):
                    chosen = [r for r in edges if subset == "all" or r[subset]]
                    ligands = ids if subset == "all" else [i for i in ids if i != 17]
                    exp = [float(kd[i]) for i in ligands]
                    pred = [scores[i] for i in ligands]
                    stats.append(
                        dict(
                            target=target,
                            pdb=pdb,
                            arm=arm,
                            subset=subset,
                            pairs=len(chosen),
                            credit=sum(r["credit"] for r in chosen),
                            spearman=(
                                spearmanr(exp, pred).statistic
                                if subset != "same_formula"
                                else ""
                            ),
                            kendall=(
                                kendalltau(exp, pred).statistic
                                if subset != "same_formula"
                                else ""
                            ),
                        )
                    )
    # The native rows must reproduce the original (independently generated) ledgers.
    truth = {
        (r["target"], r["scheme"], r["subset"]): (int(r["n"]), float(r["correct"]))
        for r in read(GAL / "results/summary.csv")
        if r["arm"] == "native"
        and r["field"] == "top"
        and r["cohort"] == "all"
        and r["scheme"] in ("0", "1", "mean")
    }
    for r in stats:
        if r["arm"] == "native":
            scheme = (
                "mean"
                if r["pdb"] == "mean"
                else str(TARGETS[r["target"]].index(r["pdb"]))
            )
            assert (r["pairs"], r["credit"]) == truth[r["target"], scheme, r["subset"]]
    bootstrap = []
    for target in TARGETS:

        def values(arm, subset):
            return {
                (r["stronger"], r["weaker"]): r["credit"]
                for r in pairs
                if r["target"] == target
                and r["pdb"] == "mean"
                and r["arm"] == arm
                and (subset == "all" or r[subset])
            }

        for subset in ("all", "triazoles"):
            compounds = ids if subset == "all" else [c for c in ids if c != 17]
            native = values("native", subset)
            left, right = np.asarray(list(native)).T
            # Resample whole compounds and preserve every receptor/seed for that compound.
            draws = np.random.default_rng(20261126).choice(
                compounds, (4000, len(compounds))
            )
            counts = np.stack(
                [(draws == i).sum(axis=1) for i in range(max(ids) + 1)], axis=1
            )
            weights = counts[:, left] * counts[:, right]
            weights = weights[weights.sum(axis=1) > 0]
            for arm in (a for a in arms if a != "native"):
                new = values(arm, subset)
                delta = np.array([new[k] - v for k, v in native.items()])
                lo, hi = np.quantile(
                    (weights @ delta) / weights.sum(axis=1), [0.025, 0.975]
                )
                bootstrap.append(
                    dict(
                        target=target,
                        subset=subset,
                        arm=arm,
                        delta_fraction=np.mean(delta),
                        low=lo,
                        high=hi,
                        valid=len(weights),
                    )
                )
    write_table(out / "bootstrap.csv", bootstrap)
    write_table(out / "geometry.csv", geometry)
    write_table(out / "affinity.csv", stats)
    write_table(out / "pairs.csv.xz", pairs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "work", type=Path, help="extracted galectin_affinity/evidence.tar.xz"
    )
    parser.add_argument("binary", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    rows = score(args.work, args.binary.resolve())
    seconds = time.perf_counter() - start
    write(args.out / "modes.csv.xz", rows)
    evaluate(rows, args.work, args.out)
    print(
        f"{len(rows)} fixed modes in {seconds:.3f}s, maximum INTER rounding discrepancy "
        f"{max(abs(r['inter_rounding_error']) for r in rows):.7f}"
    )
