"""Prepare or replay the three Aurora A ligand probes from the frozen protocol."""

import argparse
import csv
import hashlib
import itertools
import json
import resource
import subprocess
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parents[1] / "torsion_kinase_affinity/evidence.tar.xz"
RECEPTORS = ("2J4Z", "2X6E")
SEEDS = (60913, 60943, 61211, 61231)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compounds():
    with (HERE / "compounds.csv").open() as source:
        return list(csv.DictReader(source))


def prepare(work):
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem
    from rdkit.Chem import AllChem

    out = work / "prepared"
    out.mkdir(exist_ok=True, parents=True)
    if (
        digest(ARCHIVE)
        != "72c45e52eda2b79d5241b00b4efe79b2e2c19dc100f01a3ad52356d9c50cf56e"
    ):
        raise ValueError("Unexpected receptor archive")
    with tarfile.open(ARCHIVE) as archive:
        for state, ext in itertools.product(RECEPTORS, ("json", "pdbqt")):
            name = f"{state}.{ext}"
            (out / name).write_bytes(archive.extractfile(f"prepared/{name}").read())
    for row in compounds():
        mol = Chem.MolFromSmiles(row["smiles"])
        if Chem.MolToInchiKey(mol) != row["inchikey"] or Chem.GetFormalCharge(mol) != 0:
            raise ValueError(row)
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 2026
        ids = AllChem.EmbedMultipleConfs(mol, numConfs=8, params=params)
        optimized = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=1500)
        _, best = min(
            (energy, idx)
            for idx, (status, energy) in zip(ids, optimized)
            if status == 0
        )
        conf = Chem.Conformer(mol.GetConformer(best))
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        mol.SetProp("_Name", row["compound"])
        Chem.MolToMolFile(mol, str(out / (row["compound"] + ".sdf")))
        (setup,) = MoleculePreparation().prepare(mol)
        pdbqt, ok, error = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
        if not ok:
            raise ValueError(error)
        (out / (row["compound"] + ".pdbqt")).write_text(pdbqt)


def run(work, pro, official):
    if (
        digest(official)
        != "f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644"
    ):
        raise ValueError("Use the pinned official AutoDock Vina 1.2.7 binary")
    out = work / "pools"
    out.mkdir(exist_ok=True)
    names = [r["compound"] for r in compounds()]
    methods = ("pro", "official", "double")
    with (work / "runs.csv").open("w") as file:
        w = csv.writer(file, lineterminator="\n")
        w.writerow(
            (
                "method",
                "receptor",
                "compound",
                "seed",
                "cap",
                "wall_s",
                "cpu_s",
                "modes",
                "score",
                "pool_sha256",
            )
        )
        for i, (rec, name, seed) in enumerate(
            itertools.product(RECEPTORS, names, SEEDS)
        ):
            center = json.loads((work / f"prepared/{rec}.json").read_text())["center"]
            for method in methods[i % 3 :] + methods[: i % 3]:
                base = f"{method}_{rec}_{name}_{seed}"
                result = Path("pools") / (base + ".pdbqt")
                cap = 50000 if method == "double" else 25000
                args = [
                    "--receptor",
                    f"prepared/{rec}.pdbqt",
                    "--ligand",
                    f"prepared/{name}.pdbqt",
                    "--seed",
                    str(seed),
                    "--exhaustiveness",
                    "2",
                    "--cpu",
                    "1",
                    "--max_evals",
                    str(cap),
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--out",
                    str(result),
                ]
                for dim, c in zip("xyz", center):
                    args.extend([f"--center_{dim}", str(c), f"--size_{dim}", "26"])
                binary = official if method == "official" else pro
                if method != "official":
                    args.extend(["--pose_selection", "late", "--search_counters"])
                start = time.perf_counter()
                before = resource.getrusage(resource.RUSAGE_CHILDREN)
                with (out / (base + ".log")).open("w") as log:
                    subprocess.run(
                        [binary, *args],
                        cwd=work,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        check=True,
                    )
                after = resource.getrusage(resource.RUSAGE_CHILDREN)
                wall = time.perf_counter() - start
                cpu = (after.ru_utime - before.ru_utime) + (
                    after.ru_stime - before.ru_stime
                )
                scores = [
                    float(l.split()[3])
                    for l in (work / result).read_text().splitlines()
                    if l.startswith("REMARK VINA RESULT:")
                ]
                w.writerow(
                    [
                        method,
                        rec,
                        name,
                        seed,
                        cap,
                        round(wall, 6),
                        round(cpu, 6),
                        len(scores),
                        scores[0],
                        digest(work / result),
                    ]
                )
                file.flush()
    (work / "binaries.json").write_text(
        json.dumps({"pro": digest(pro), "official": digest(official)}, indent=2) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True, dest="stage")
    p = sub.add_parser("prepare")
    p.add_argument("work", type=Path)
    p = sub.add_parser("run")
    p.add_argument("work", type=Path)
    p.add_argument("pro", type=Path)
    p.add_argument("official", type=Path)
    a = parser.parse_args()
    if a.stage == "prepare":
        prepare(a.work.resolve())
    else:
        run(a.work.resolve(), a.pro.resolve(), a.official.resolve())
