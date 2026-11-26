"""Rebuild four receptor states and thirteen endpoint-independent ligand inputs."""

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from urllib.request import urlopen

os.environ.setdefault("OPENMM_CPU_THREADS", "1")
import numpy as np
from openmm import Platform
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem

# Reuse the branch's current deterministic ETKDG/MMFF/Meeko ligand preparation.
sys.path.insert(
    0, str(Path(__file__).resolve().parents[1] / "kinase_state_affinity" / "source")
)
from prepare import ligand

STATES = {
    "1GZW": ("Gal1", "B", "D", ("GAL", "BGC")),
    "4Q26": ("Gal1", "A", "C", ("GAL", "NDG")),
    "1KJL": ("Gal3", "A", "B", ("GAL", "NAG")),
    "3T1L": ("Gal3", "A", "A", ("MQT",)),
}
ROOT = Path(__file__).resolve().parent
GFA = "https://raw.githubusercontent.com/choutkaj/GFA-ligand-set/04dd9b0024ba3a81730731d256b6e7b285977aec/GFA.csv"
PAPER = "10.1021/acs.jmedchem.7b01626"


def download(path, url):
    if not path.exists():
        path.write_bytes(urlopen(url, timeout=60).read())
    expected = json.loads((ROOT / "raw_sha256.json").read_text())[path.name]
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Changed source: {path.name}")


def prepare_receptor(work, pdb):
    target, chain, donor_chain, components = STATES[pdb]
    prepared = work / "prepared"
    rows = (work / "raw" / f"{pdb}.pdb").read_text().splitlines()
    donor = []
    protein = []
    for line in rows:
        if (
            line[:6] not in ("HETATM", "ATOM  ")
            or line[16] not in " A"
            or line[76:78].strip() in ("H", "D")
        ):
            continue
        comp, atom = line[17:20].strip(), line[12:16].strip()
        if line.startswith("HETATM") and comp in components and line[21] == donor_chain:
            donor.append([float(line[k : k + 8]) for k in (30, 38, 46)])
        if line[21] != chain:
            continue
        if (
            line.startswith("HETATM")
            and comp == "CME"
            and atom in ("N", "CA", "C", "O", "CB", "SG")
        ):
            line = "ATOM  " + line[6:17] + "CYS" + line[20:]
        if line.startswith("ATOM"):
            protein.append(line[:16] + " " + line[17:])
    if len(donor) < 12 or not protein:
        raise ValueError(pdb)
    center = np.mean(donor, axis=0)
    (prepared / f"{pdb}.box").write_text(
        "".join(
            f"center_{axis} = {v:.5f}\nsize_{axis} = 24\n"
            for axis, v in zip("xyz", center)
        )
    )
    p = prepared / f"{pdb}-protein.pdb"
    p.write_text("\n".join(protein) + "\nEND\n")
    fixer = PDBFixer(filename=str(p), platform=Platform.getPlatformByName("CPU"))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    missing = {
        f"{r.name}:{r.id}": [a.name for a in atoms]
        for r, atoms in fixer.missingAtoms.items()
    }
    (prepared / f"{pdb}-missing.json").write_text(json.dumps(missing, indent=2) + "\n")
    fixer.addMissingAtoms(seed=2026)
    with p.open("w") as out:
        PDBFile.writeFile(fixer.topology, fixer.positions, out, keepIds=True)
    args = [
        sys.executable,
        "-m",
        "meeko.cli.mk_prepare_receptor",
        "--read_pdb",
        p,
        "-p",
        prepared / f"{pdb}.pdbqt",
    ]
    if target == "Gal3":
        args.extend(["--set_template", "A:158=HID"])
    with (prepared / f"{pdb}.log").open("w") as log:
        subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, check=True)


def prepare(work):
    for name in ("raw", "prepared"):
        (work / name).mkdir(parents=True, exist_ok=True)
    source = work / "raw" / "GFA.csv"
    download(source, GFA)
    records = {
        int(r["Name in source"]): r
        for r in csv.DictReader(source.open(encoding="utf-8-sig"))
        if r["Citation"] == PAPER
        and r["Name in source"].isdigit()
        and 5 <= int(r["Name in source"]) <= 17
    }
    reference = list(csv.DictReader((ROOT / "table1.csv").open()))
    with (work / "panel.csv").open("w") as outf:
        writer = csv.writer(outf)
        writer.writerow(
            [
                "compound",
                "gfa_id",
                "Gal1_Kd_uM",
                "Gal3_Kd_uM",
                "heavy_atoms",
                "fluorines",
                "smiles",
            ]
        )
        for row in reference:
            cid = row["compound"]
            rec = records[int(cid)]
            for number in (1, 3):
                if rec[f"Gal-{number} Kd operator"] != "=" or float(
                    rec[f"Gal-{number} Kd value (μM)"]
                ) != float(row[f"Gal{number}_Kd_uM"]):
                    raise ValueError(f"Endpoint mismatch: {cid}")
            smiles = rec["SMILES"]
            mol = Chem.MolFromSmiles(smiles)
            writer.writerow(
                [
                    cid,
                    rec["GFA_ID"],
                    row["Gal1_Kd_uM"],
                    row["Gal3_Kd_uM"],
                    mol.GetNumHeavyAtoms(),
                    sum(a.GetAtomicNum() == 9 for a in mol.GetAtoms()),
                    smiles,
                ]
            )
            ligand(work, {"smiles": smiles, "cid": cid})
    for pdb in STATES:
        download(
            work / "raw" / f"{pdb}.pdb", f"https://files.rcsb.org/download/{pdb}.pdb"
        )
        prepare_receptor(work, pdb)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    prepare(parser.parse_args().work)
