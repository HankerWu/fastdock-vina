"""Replay a same-assay Factor Xa Ki series using an independent, common receptor."""

import argparse
import csv
import hashlib
import itertools
import json
import subprocess
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from statistics import mean

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent
ACTIVITY = "https://www.ebi.ac.uk/chembl/api/data/activity.json?assay_chembl_id=CHEMBL661804&limit=1000"
RECEPTORS = {
    "2P16": (
        "GG2",
        34,
        "00ee06bbd8635d732edc7a5acd3ba255d8489976a8a382fd51cf295cdb229905",
    ),
    "1FJS": (
        "Z34",
        38,
        "240307991202f732540d06fe585a8291d541412d410a2539d5d8ce6151a6d53e",
    ),
}
SEEDS = (271828, 314159)


def write(path, rows):
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read(path):
    with path.open() as inp:
        return list(csv.DictReader(inp))


def curate(api):
    data = json.loads(api)["activities"]
    rows = []
    for r in sorted(data, key=lambda row: row["activity_id"]):
        if (
            r["standard_relation"],
            r["standard_type"],
            r["standard_units"],
            r["document_chembl_id"],
            r["assay_chembl_id"],
            r["target_chembl_id"],
        ) != ("=", "Ki", "nM", "CHEMBL1133358", "CHEMBL661804", "CHEMBL244"):
            raise ValueError("Assay changed, check source")
        mol = Chem.MolFromSmiles(r["canonical_smiles"])
        unresolved = any(
            tag == "?"
            for _, tag in Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        )
        # Fixed dominant amidinium and carboxylate states for this one series.
        # Both the primary aryl amidine and N-acetimidoylpiperidine are strong bases.
        amidines = mol.GetSubstructMatches(Chem.MolFromSmarts("[C;X3](=[NH1])-[N;X3]"))
        acids = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[OH1]"))
        if len(amidines) != 2 or len(acids) > 1:
            raise ValueError(f"Unexpected chemistry in {r['molecule_chembl_id']}")
        for _, nitrogen, _ in amidines:
            n = mol.GetAtomWithIdx(nitrogen)
            n.SetFormalCharge(1)
            n.SetNumExplicitHs(2)
            n.SetNoImplicit(True)
        for _, _, oxygen in acids:
            o = mol.GetAtomWithIdx(oxygen)
            o.SetFormalCharge(-1)
            o.SetNumExplicitHs(0)
            o.SetNoImplicit(True)
        Chem.SanitizeMol(mol)
        rows.append(
            {
                "ligand": r["molecule_chembl_id"],
                "activity": r["activity_id"],
                "ki_nm": r["standard_value"],
                "original": r["canonical_smiles"],
                "docked_smiles": Chem.MolToSmiles(mol),
                "heavy": mol.GetNumHeavyAtoms(),
                "charge": Chem.GetFormalCharge(mol),
                "included": int(not unresolved),
            }
        )
    assert len({r["ligand"] for r in rows}) == len(rows)
    return rows


def prepare(work, frozen, rec):
    component, expected_atoms, sha256 = RECEPTORS[rec]
    source = work / f"{rec}.pdb"
    if not source.exists():
        url = f"https://files.rcsb.org/download/{rec}.pdb"
        source.write_bytes(urllib.request.urlopen(url, timeout=60).read())
    if hashlib.sha256(source.read_bytes()).hexdigest() != sha256:
        raise ValueError(f"{rec} source digest changed; verify structure")
    lines = source.read_text().splitlines()
    atoms = [
        l[:16] + " " + l[17:] for l in lines if l.startswith("ATOM  ") and l[16] in " A"
    ]
    protein = work / "protein.pdb"
    protein.write_text("\n".join(atoms) + "\nEND\n")
    with (work / "prep.log").open("w") as log:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "meeko.cli.mk_prepare_receptor",
                "--read_pdb",
                protein,
                "-p",
                work / f"{rec}.pdbqt",
            ],
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    bound_ligand = [
        l
        for l in lines
        if l.startswith("HETATM") and l[17:20] == component and l[21] == "A"
    ]
    assert len(bound_ligand) == expected_atoms
    center = [mean(float(l[t : t + 8]) for l in bound_ligand) for t in (30, 38, 46)]
    (work / "box.txt").write_text(
        "".join(
            f"center_{axis} = {v:.4f}\nsize_{axis} = 26\n"
            for axis, v in zip("xyz", center)
        )
    )
    for r in frozen:
        if int(r["included"]):
            m = Chem.AddHs(Chem.MolFromSmiles(r["docked_smiles"]))
            options = AllChem.ETKDGv3()
            options.randomSeed = 2026
            if AllChem.EmbedMolecule(m, options) != 0:
                raise ValueError(r["ligand"])
            # A common, label-independent preparation for all ligands.
            status = AllChem.MMFFOptimizeMolecule(m, maxIters=1000)
            if status != 0:
                raise ValueError(f"MMFF did not converge: {r['ligand']}")
            (setup,) = MoleculePreparation().prepare(m)
            pdbqt, ok, reason = PDBQTWriterLegacy.write_string(
                setup, add_index_map=True
            )
            if not ok:
                raise ValueError(reason)
            (work / f"{r['ligand']}.pdbqt").write_text(pdbqt)


def dock(work, binary, jobs, frozen, rec):
    def one(case):
        r, seed = case
        base = work / f"{r['ligand']}_{seed}"
        with base.with_suffix(".log").open("w") as log:
            subprocess.run(
                [
                    binary,
                    "--receptor",
                    work / f"{rec}.pdbqt",
                    "--ligand",
                    work / f"{r['ligand']}.pdbqt",
                    "--config",
                    work / "box.txt",
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
                    "--out",
                    base.with_suffix(".pdbqt"),
                ],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
            )

    cases = itertools.product([r for r in frozen if int(r["included"])], SEEDS)
    with ThreadPoolExecutor(jobs) as pool:
        list(pool.map(one, cases))


def reported_scores(lines):
    """Use the scores at deposited output precision on one common model."""

    def energy(prefix):
        values = [
            float(line.removeprefix(prefix).split()[0])
            for line in lines
            if line.startswith(prefix)
        ]
        (value,) = values
        return value

    return {
        "native": energy("REMARK VINA RESULT:"),
        "no_divisor": round(
            energy("REMARK INTER + INTRA:") - energy("REMARK UNBOUND:"), 3
        ),
    }


def collect(work, frozen, output, rec):
    out = []
    for r in frozen:
        if not int(r["included"]):
            continue
        for seed in SEEDS:
            p = work / f"{r['ligand']}_{seed}.pdbqt"
            lines = p.read_text().split("ENDMDL", 1)[0].splitlines()

            out.append(
                dict(
                    target="fxa",
                    ligand=r["ligand"],
                    receptor=rec,
                    seed=seed,
                    measurement=r["ki_nm"],
                    heavy=r["heavy"],
                    charge=r["charge"],
                    pose_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                    **reported_scores(lines),
                )
            )
    write(output, out)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("fetch", "prepare", "dock", "collect"))
    p.add_argument("work", type=Path)
    p.add_argument("--binary", type=Path)
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--receptor", choices=RECEPTORS, default="2P16")
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    table = HERE / "fxa_series.csv"
    if a.action == "fetch":
        write(table, curate(urllib.request.urlopen(ACTIVITY, timeout=60).read()))
    elif a.action == "prepare":
        prepare(a.work, read(table), a.receptor)
    elif a.action == "dock":
        dock(a.work, a.binary.resolve(), a.jobs, read(table), a.receptor)
    else:
        collect(
            a.work,
            read(table),
            HERE / "results" / f"fxa_{a.receptor}_cells.csv",
            a.receptor,
        )


if __name__ == "__main__":
    main()
