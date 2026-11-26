"""Hold the optimizer fixed while comparing BPTF water models in measured Kd assays."""
import argparse
import csv
import itertools
import json
import statistics
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
ARMS = ("dry", "steric", "polar")
SEEDS = {48000: (84551, 84559, 84571), 120000: (93701, 93719, 93731)}


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def prepare(work):
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    sys.path.insert(0, str(ROOT / "benchmarks/kinase_state_affinity/source"))
    from prepare import ligand
    definitions = json.loads((HERE / "identities.json").read_text())
    rejected = {r["key"] for r in definitions["donors_and_receivers"]}
    (work / "prepared").mkdir(parents=True, exist_ok=True)
    for row in definitions["ligands"]:
        mol = Chem.MolFromSmiles(row["smiles"])
        key = Chem.MolToInchiKey(mol)
        assert key == row["key"] and key not in rejected
        assert rdMolDescriptors.CalcMolFormula(mol) == row["neutral_formula"]
        if row["protonate"]:
            # Primary/tertiary aliphatic amine, distinct from conjugated anilines/
            # sulfonamides and from the nitrogens in the aromatic ring systems.
            pattern = Chem.MolFromSmarts("[N;X3;+0;!$(N-[a]);!$(N-C=[O,N,S]);!$(N-S)]")
            (i,), = mol.GetSubstructMatches(pattern)
            atom = mol.GetAtomWithIdx(i)
            atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
            atom.SetNoImplicit(True)
            atom.SetFormalCharge(1)
            Chem.SanitizeMol(mol)
        ligand(work, dict(cid=row["compound"], smiles=Chem.MolToSmiles(mol)))


def cases():
    return json.loads((HERE / "identities.json").read_text())["ligands"]


def run(work, receptors, binary, jobs, arms, seeds):
    params = json.loads((HERE / "centers.json").read_text())
    pools = work / "pools"
    pools.mkdir(parents=True, exist_ok=True)

    def dock(cell):
        row, pdb, cap, index = cell
        seed = seeds[cap][index]
        for arm in arms[index:] + arms[:index]:
            stem = pools / f"{row['compound']}_{pdb}_{cap}_{seed}_{arm}"
            cmd = [binary, "--ligand", work / "prepared" / (row["compound"] + ".pdbqt"),
                   "--receptor", receptors / arm / pdb / "receptor.pdbqt",
                   "--pose_selection", "late", "--refinement_optimizer", "geometry", "--cpu", "1",
                   "--search_counters", "--exhaustiveness", "3", "--max_evals", str(cap),
                   "--num_modes", "9", "--energy_range", "100", "--seed", str(seed),
                   "--out", stem.with_suffix(".pdbqt")]
            for axis, center in zip("xyz", params[pdb]):
                cmd += ["--center_" + axis, str(center), "--size_" + axis, "24"]
            with stem.with_suffix(".log").open("w") as f:
                subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)

    tasks = itertools.product(cases(), params, seeds, range(3))
    with ThreadPoolExecutor(jobs) as ex:
        list(ex.map(dock, tasks))


def analyze(work, arms, seeds):
    rows = cases()
    scores, cells, pairs = {}, [], []
    centers = json.loads((HERE / "centers.json").read_text())
    for row, pdb, cap, arm in itertools.product(rows, centers, seeds, arms):
        values = []
        for seed in seeds[cap]:
            prefix = work / "pools" / f"{row['compound']}_{pdb}_{cap}_{seed}_{arm}"
            modes = [float(l.split()[3]) for l in prefix.with_suffix(".pdbqt").read_text().splitlines()
                     if l.startswith("REMARK VINA RESULT:")]
            assert modes and modes[0] == min(modes)
            values.append(modes[0])
            cells.append(dict(compound=row["compound"], receptor=pdb, cap=cap, seed=seed,
                              arm=arm, score=modes[0], modes=len(modes)))
        scores[pdb, cap, arm, row["compound"]] = statistics.median(values)
    for a, b in itertools.combinations(rows, 2):
        if a["assay"] != b["assay"] or a["kd_nM"] == b["kd_nM"]:
            continue
        a, b = sorted((a, b), key=lambda r: r["kd_nM"])
        for pdb, cap in itertools.product(centers, seeds):
            d = dict(assay=a["assay"], receptor=pdb, cap=cap, stronger=a["compound"], weaker=b["compound"],
                     fold=b["kd_nM"] / a["kd_nM"])
            for arm in arms:
                d[arm] = round(scores[pdb, cap, arm, b["compound"]] - scores[pdb, cap, arm, a["compound"]], 6)
            pairs.append(d)
    write(work / "cells.csv", cells)
    write(work / "pairs.csv", pairs)
    summaries = []
    for assay, pdb, cap in itertools.product(sorted({r["assay"] for r in rows}), centers, seeds):
        for scope in ("all", "without-24") if assay.startswith("Zahid") else ("all",):
            selected = [p for p in pairs if p["receptor"] == pdb and p["cap"] == cap and p["assay"] == assay
                        and p["fold"] >= 3 and (scope == "all" or "24" not in (p["stronger"], p["weaker"]))]
            summaries.append(dict(assay=assay, receptor=pdb, cap=cap, scope=scope, pairs=len(selected),
                                  **{arm: sum((p[arm] > 0) + 0.5 * (p[arm] == 0) for p in selected) for arm in arms}))
    write(work / "summary.csv", summaries)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--receptors", type=Path)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--confirmation", action="store_true")
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    work = args.work.resolve()
    arms, seeds = (("dry", "polar", "supported", "registered"), {120000: (95701, 95713, 95717)}) if args.confirmation else (ARMS, SEEDS)
    if args.prepare:
        prepare(work)
    if args.binary:
        run(work, args.receptors.resolve(), args.binary.resolve(), args.jobs, arms, seeds)
    analyze(work, arms, seeds)
