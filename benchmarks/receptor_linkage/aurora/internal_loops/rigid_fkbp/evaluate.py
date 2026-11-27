"""Reconstruct all pose and ranking endpoints directly from the deposited coordinates and printed pools."""

import argparse
import csv
import itertools
import json
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial.transform import Rotation
from scipy.stats import kendalltau, spearmanr

from run_panel import ARMS, BENCH, ROOT, RECEIVERS, SOURCE, prepare, sha, table

sys.path.insert(0, str(SOURCE))
from analyze import calculate_components, read_map

sys.path.insert(0, str(BENCH / "human_panel"))
from verify import reference

sys.path.insert(0, str(ROOT / "benchmarks/proposal_sweep/fkbp"))
from audit import molecule, stereo, xyz


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def crystal_positions(source, receiver, query):
    a, ca, mol = source
    b, cb, _ = receiver
    keys = sorted(
        k
        for k in a.keys() & b.keys()
        if np.linalg.norm(a[k] - ca) < 14 and np.linalg.norm(b[k] - cb) < 14
    )
    x, y = (np.array([c[k] for k in keys]) for c in (a, b))
    rot, _ = Rotation.align_vectors(y - y.mean(0), x - x.mean(0))
    points = rot.apply(mol.GetConformer().GetPositions() - x.mean(0)) + y.mean(0)
    mappings = mol.GetSubstructMatches(
        query, uniquify=False, useChirality=True, maxMatches=10001
    )
    assert mol.GetNumAtoms() == query.GetNumAtoms() and 0 < len(mappings) < 10001
    return points[np.asarray(mappings)]


def inspect(work, out, arms=ARMS, comparison="pivot"):
    out.mkdir(parents=True, exist_ok=True)
    hashes = prepare(work)
    dep = BENCH / "human_panel/depositions.tar.xz"
    assay = BENCH / "human_panel/affinity/ligands.csv"
    for p in (dep, assay, SOURCE / "results/metadata.json"):
        hashes[str(p.relative_to(ROOT))] = sha(p)
    (out / "inputs.json").write_text(json.dumps(hashes, indent=2) + "\n")
    with tarfile.open(dep) as tar:
        for name in ("1FKG.pdb", "1FKH.pdb", "SB3.cif", "SBX.cif"):
            tar.extract("raw/" + name, work, filter="data")
    refs = {
        r: reference(work, r, c) for r, c in zip(RECEIVERS, ("SB3", "SBX"), strict=True)
    }
    expected_smiles = {
        int(k["compound"]): Chem.MolToSmiles(Chem.MolFromSmiles(k["smiles"]))
        for k in read(assay)
    }
    labels = {
        int(k["compound"]): (float(k["ki_low_nm"]), float(k["ki_high_nm"]))
        for k in read(assay)
    }
    # Independent read from Table II of Holt et al., JACS 1993, p.9927.
    assert [labels[n] for n in range(4, 12)] == [
        (186, 186),
        (110, 110),
        (12, 12),
        (250, 250),
        (10, 10),
        (7, 7),
        (300, 600),
        (300, 300),
    ]
    maps = {
        r: {
            t: read_map(work / "maps" / r / f"receptor.{t}.map")
            for t in ("A", "C", "N", "OA", "e", "d")
        }
        for r in RECEIVERS
    }
    runs = [row for arm in arms for row in read(work / f"runs-{arm}.csv")]
    replay = read(work / "counter_replay.csv")
    by_cell = {(r["receiver"], r["compound"]): r for r in replay}
    assert len(by_cell) == 16
    table(out / "counter_replay.csv", replay)
    for row in runs:
        r, n, arm = row["receiver"], row["compound"], row["arm"]
        stem = f"{r}_{n}"
        if arm == "official":
            counterpart = by_cell[r, n]
            counter_dir = work / "official_counter"
            equal = (counter_dir / (stem + ".pdbqt")).read_bytes() == (
                work / "pools/official" / (stem + ".pdbqt")
            ).read_bytes()
            assert equal == bool(int(counterpart["byte_equal"]))
            if not equal:
                continue
        else:
            counter_dir = work / "pools" / arm
            counterpart = row
        log = (counter_dir / (stem + ".log")).read_text()
        for field, tag in (
            ("hunt_evals", "hunt_evals"),
            ("second_evals", "second_evals"),
            ("steps", "hunt_calls"),
        ):
            terms = list(map(int, re.findall(r"\b" + tag + r"=(\d+)", log)))
            assert len(terms) == 3 and sum(terms) == int(counterpart[field])
            row[field] = sum(terms)
    assert (
        len({(r["arm"], r["receiver"], r["compound"]) for r in runs}) == len(runs) == len(arms) * 16
    )
    modes, poses, desc, scores, costs = [], [], {}, {}, []
    rescored = {}
    for run in runs:
        arm, rec, n = run["arm"], run["receiver"], int(run["compound"])
        ligand = (work / "affinity" / f"{n}.pdbqt").read_text()
        mol, serial, bonds = molecule(ligand)
        typed = Chem.MolToSmiles(mol)
        assert typed == expected_smiles[n]
        inp = xyz(ligand, serial)
        assert stereo(mol, inp) == typed
        desc[n] = dict(
            compound=n,
            heavy_atoms=len(serial),
            formal_charge=Chem.GetFormalCharge(mol),
            torsdof=int(
                next(
                    l.split()[1] for l in ligand.splitlines() if l.startswith("TORSDOF")
                )
            ),
        )

        def lengths(x):
            return np.linalg.norm(x[bonds[:, 0]] - x[bonds[:, 1]], axis=1)

        initial_bonds = lengths(inp)
        crystal = {8: "1FKG", 9: "1FKH"}.get(n)
        if crystal:
            fitted = crystal_positions(refs[crystal], refs[rec], mol)
        p = work / "pools" / arm / f"{rec}_{n}.pdbqt"
        assert sha(p) == run["sha256"]
        log = p.with_suffix(".log").read_text()
        assert "(random seed: 524287)" in log and "Scoring function : ad4" in log
        blocks = [b for b in p.read_text().split("ENDMDL") if b.strip()]
        assert len(blocks) == int(run["modes"]) and len(blocks) > 0
        rmsds = []
        for rank, block in enumerate(blocks, 1):
            lines = block.splitlines()
            output_mol, output_serial, output_bonds = molecule(block)
            assert Chem.MolToSmiles(output_mol) == typed and output_serial == serial
            assert np.array_equal(output_bonds, bonds)
            positions = xyz(block, serial)
            reported = float(
                next(l for l in lines if l.startswith("REMARK VINA RESULT:")).split()[3]
            )
            operator = "separate" if arm in ("official", "official_high") else "combined"
            rebuilt = sum(calculate_components(lines, maps[rec], operator=operator).values())
            changed = stereo(mol, positions) != typed
            drift = float(np.max(np.abs(lengths(positions) - initial_bonds)))
            if rank == 1:
                scores[arm, rec, n] = reported
                rescored[arm, rec, n] = rebuilt
            error = ""
            if crystal:
                error = float(
                    np.min(
                        np.sqrt(
                            np.mean(np.sum((fitted - positions) ** 2, axis=2), axis=1)
                        )
                    )
                )
                rmsds.append(error)
            modes.append(
                dict(
                    arm=arm,
                    receiver=rec,
                    compound=n,
                    rank=rank,
                    reported=reported,
                    ad4_from_xyz=rebuilt,
                    rmsd_A=error,
                    stereo_changed=int(changed),
                    max_bond_delta_A=drift,
                )
            )
        if crystal:
            poses.append(
                dict(
                    arm=arm,
                    receiver=rec,
                    compound=n,
                    crystal=crystal,
                    top1=rmsds[0],
                    best5=min(rmsds[:5]),
                    oracle=min(rmsds),
                    regret=rmsds[0] - min(rmsds),
                )
            )
    table(out / "runs.csv", runs)
    table(out / "modes.csv", modes)
    table(out / "pose_cells.csv", poses)
    assays, pairs, individual, distributions = [], [], [], []
    for arm in arms:
        picked = [r for r in poses if r["arm"] == arm]
        values = [r["top1"] for r in picked]
        qs = np.quantile(values, [0.25, 0.5, 0.75, 0.9, 0.95])
        distributions.append(
            dict(
                arm=arm,
                cells=len(values),
                mean=np.mean(values),
                **dict(zip(("q25", "median", "q75", "p90", "p95"), qs, strict=True)),
                **{
                    k: np.mean([r[k] for r in picked])
                    for k in ("best5", "oracle", "regret")
                },
                gt5=sum(v > 5 for v in values),
                ge8=sum(v >= 8 for v in values),
                gt10=sum(v > 10 for v in values),
            )
        )
        group = [r for r in runs if r["arm"] == arm]
        details = [m for m in modes if m["arm"] == arm]
        costs.append(
            dict(
                arm=arm,
                runs=len(group),
                empty=sum(int(r["modes"]) == 0 for r in group),
                modes=len(details),
                stereo_changed=sum(r["stereo_changed"] for r in details),
                worst_bond_delta_A=max(r["max_bond_delta_A"] for r in details),
                worst_score_xyz_delta=max(
                    abs(r["reported"] - r["ad4_from_xyz"]) for r in details
                ),
                wall_s=sum(float(r["wall_s"]) for r in group),
                cpu_s=sum(float(r["cpu_s"]) for r in group),
                actual_evals=(
                    sum(int(r["hunt_evals"]) + int(r["second_evals"]) for r in group)
                    if all(r["hunt_evals"] != "" for r in group)
                    else ""
                ),
            )
        )
        for n in labels:
            scores[arm, "mean", n] = float(
                np.mean([scores[arm, r, n] for r in RECEIVERS])
            )
            rescored[arm, "mean", n] = float(
                np.mean([rescored[arm, r, n] for r in RECEIVERS])
            )
        for rec in (*RECEIVERS, "mean"):
            subset = {n: scores[arm, rec, n] for n in labels}
            for n, sc in subset.items():
                individual.append(
                    dict(
                        arm=arm,
                        receiver=rec,
                        compound=n,
                        top1_score=sc,
                        ki_low_nM=labels[n][0],
                        ki_high_nM=labels[n][1],
                    )
                )
            credits = []
            for i, j in itertools.permutations(labels, 2):
                if labels[j][0] < 3 * labels[i][1]:
                    continue
                margin = subset[j] - subset[i]
                credit = 1 if margin > 1e-9 else 0 if margin < -1e-9 else 0.5
                credits.append(credit)
                pairs.append(
                    dict(
                        arm=arm,
                        receiver=rec,
                        stronger=i,
                        weaker=j,
                        minimum_fold=labels[j][0] / labels[i][1],
                        predicted_margin=margin,
                        credit=credit,
                        xyz_margin=rescored[arm, rec, j] - rescored[arm, rec, i],
                    )
                )
            assert len(credits) == 15
            points = [n for n, (lo, hi) in labels.items() if lo == hi]
            reported = [labels[n][0] for n in points]
            pred = [subset[n] for n in points]
            assays.append(
                dict(
                    arm=arm,
                    receiver=rec,
                    correct=sum(credits),
                    denominator=len(credits),
                    points=len(points),
                    spearman=spearmanr(reported, pred).statistic,
                    kendall_b=kendalltau(reported, pred).statistic,
                    cliff8_10=subset[10] - subset[8],
                    cliff9_11=subset[11] - subset[9],
                )
            )
    contrast = f"score_change_{comparison}_minus_base"
    for n in labels:
        desc[n][contrast] = (
            scores[comparison, "mean", n] - scores["base", "mean", n]
        )
    for name, data in (
        ("assay", assays),
        ("pairs", pairs),
        ("scores", individual),
        ("pose_distributions", distributions),
        ("cost_chemistry", costs),
        ("descriptors", list(desc.values())),
    ):
        table(out / f"{name}.csv", data)
    bias = {
        k: spearmanr(
            [d[k] for d in desc.values()],
            [d[contrast] for d in desc.values()],
        ).statistic
        for k in ("heavy_atoms", "torsdof")
    }
    assert max(r["worst_score_xyz_delta"] for r in costs) < 0.03
    assert all(r["credit"] == float(r["xyz_margin"] > 0) for r in pairs)
    (out / "bias.json").write_text(json.dumps(bias, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    inspect(args.work.resolve(), args.out.resolve())
