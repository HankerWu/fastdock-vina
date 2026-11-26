"""Holt et al. (1993), Table II. Compare a priori >=3-fold Ki,app orders."""

import argparse
import csv
import itertools
import json
import statistics
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from meeko import MoleculePreparation, PDBQTWriterLegacy
from panel import CAP, SEEDS
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.stats import spearmanr

# The ester alcohol changes; shared acid 30 has the verified (2S) center.
ACID = "CCC(C)(C)C(=O)C(=O)N1CCCC[C@H]1C(=O)O"
SERIES = [
    (4, "CCCC1CCCCC1", 186, 186, 24),
    (5, "CCCc1ccccc1", 110, 110, 20),
    (6, "CCCc1cc(OC)c(OC)c(OC)c1", 12, 12, 5),
    (7, "[C@H](CCc1ccccc1)C(C)(C)C=C", 250, 250, 30),
    (8, "[C@H](CCc1ccccc1)c1ccccc1", 10, 10, 1),
    (9, "[C@H](CCc1ccccc1)C1CCCCC1", 7, 7, 2),
    (10, "[C@@H](CCc1ccccc1)c1ccccc1", 300, 600, None),
    (11, "[C@@H](CCc1ccccc1)C1CCCCC1", 300, 300, 40),
]


def table(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def prepare(work, aff):
    labels = []
    for compound, sub, lo, hi, err in SERIES:
        mol = Chem.MolFromSmiles(ACID + sub)
        smi = Chem.MolToSmiles(mol)
        mol = Chem.MolFromSmiles(smi)
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        assert sum(c == "S" for _, c in centers) == (2 if compound in (10, 11) else 1)
        if compound in (7, 8, 9):
            assert sum(c == "R" for _, c in centers) == 1
        if compound in (8, 9):
            pdb = {8: "1FKG", 9: "1FKH"}[compound]
            assert (
                smi
                == json.loads((work / pdb / "preparation.json").read_text())[
                    "canonical_smiles"
                ]
            )
        m = Chem.AddHs(mol)
        p = AllChem.ETKDGv3()
        p.randomSeed = 43103
        assert (
            AllChem.EmbedMolecule(m, p) == 0
            and AllChem.MMFFOptimizeMolecule(m, maxIters=1000) == 0
        )
        (setup,) = MoleculePreparation().prepare(m)
        qt, success, msg = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
        if not success:
            raise ValueError(msg)
        (aff / f"{compound}.pdbqt").write_text(qt)
        labels.append(
            {
                "compound": compound,
                "smiles": smi,
                "heavy_atoms": mol.GetNumHeavyAtoms(),
                "ki_low_nm": lo,
                "ki_high_nm": hi,
                "reported_error_nm": err,
            }
        )
    table(aff / "ligands.csv", labels)
    return labels


def run(work, aff, binary, jobs, seeds):
    labels = prepare(work, aff)

    def one(task):
        label, receiver, seed, arm = task
        n = label["compound"]
        name = f"{n}_{receiver}_{seed}_{arm}"
        pool = aff / "runs" / f"{name}.pdbqt"
        log = pool.with_suffix(".log")
        pool.parent.mkdir(exist_ok=True)
        cmd = [
            str(binary),
            "--receptor",
            str(work / receiver / "receptor.pdbqt"),
            "--ligand",
            str(aff / f"{n}.pdbqt"),
            "--pose_selection",
            "late",
            "--cpu",
            "1",
            "--exhaustiveness",
            "3",
            "--max_evals",
            str(CAP),
            "--num_modes",
            "9",
            "--refine_candidates",
            "15",
            "--energy_range",
            "100",
            "--search_counters",
            "--seed",
            str(seed),
            "--out",
            str(pool),
        ]
        for axis, v in zip(
            "xyz",
            json.loads((work / receiver / "preparation.json").read_text())["center"],
        ):
            cmd.extend(["--center_" + axis, str(v), "--size_" + axis, "24"])
        if arm == "cubic":
            cmd.append("--cubic_search_grid")
        with log.open("w") as f:
            subprocess.run(
                cmd, stdout=f, stderr=subprocess.STDOUT, check=True, timeout=300
            )
        score = float(
            next(
                l
                for l in pool.read_text().splitlines()
                if l.startswith("REMARK VINA RESULT:")
            ).split()[3]
        )
        return {
            "compound": n,
            "receiver": receiver,
            "seed": seed,
            "arm": arm,
            "score": score,
        }

    with ThreadPoolExecutor(jobs) as pool:
        cells = list(
            pool.map(
                one,
                itertools.product(labels, ["1FKG", "1FKH"], seeds, ["native", "cubic"]),
            )
        )
    table(aff / "cells.csv", cells)
    scores = {}
    summ = []
    for label, arm in itertools.product(labels, ["native", "cubic"]):
        n = label["compound"]
        meds = []
        for r in ["1FKG", "1FKH"]:
            meds.append(
                statistics.median(
                    row["score"]
                    for row in cells
                    if row["arm"] == arm
                    and row["receiver"] == r
                    and row["compound"] == n
                )
            )
        scores[n, arm] = statistics.mean(meds)
        summ.append(
            {
                "compound": n,
                "arm": arm,
                "median_1fkg": meds[0],
                "median_1fkh": meds[1],
                "mean": scores[n, arm],
            }
        )
    table(aff / "scores.csv", summ)
    pairs = []
    for x, y in itertools.permutations(labels, 2):
        if y["ki_low_nm"] < 3 * x["ki_high_nm"]:
            continue
        i, j = x["compound"], y["compound"]
        row = {
            "stronger": i,
            "weaker": j,
            "minimum_fold": y["ki_low_nm"] / x["ki_high_nm"],
        }
        for arm in ["native", "cubic"]:
            diff = scores[j, arm] - scores[i, arm]
            row[arm + "_delta"] = diff
            row[arm + "_credit"] = (
                1 if diff > 1e-10 else 0.5 if abs(diff) < 1e-10 else 0
            )
        pairs.append(row)
    table(aff / "pairs.csv", pairs)
    corr = []
    for arm in ["native", "cubic"]:
        sels = [x for x in labels if x["ki_low_nm"] == x["ki_high_nm"]]
        rho = spearmanr(
            [x["ki_low_nm"] for x in sels], [scores[x["compound"], arm] for x in sels]
        ).statistic
        corr.append(
            {
                "arm": arm,
                "ki_pair_count": len(pairs),
                "pair_credit": sum(x[arm + "_credit"] for x in pairs),
                "rho_seven_point_ki": rho,
            }
        )
    table(aff / "summary.csv", corr)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--confirmation", action="store_true")
    a = p.parse_args()
    work = a.work.resolve()
    aff = work / ("affinity_confirmation" if a.confirmation else "affinity")
    aff.mkdir(exist_ok=True)
    run(
        work,
        aff,
        a.binary.resolve(),
        a.jobs,
        [92461, 148721, 211693] if a.confirmation else SEEDS,
    )
