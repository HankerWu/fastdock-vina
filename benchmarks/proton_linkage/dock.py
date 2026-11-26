"""A paired native/point-anchor dock on bovine CA II and seven independently titrated ligands."""

import argparse
import csv
import hashlib
import itertools
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "benchmarks/anchor_search"))
from holdout import receptor

HASHES = {
    "1V9E": "4677d09c900e324c7615cd0300865a0bc6cca3e3979757d94ebeed117086f591",
    "6SKT": "3c9886882f2081ab43070de90493d49726cacc054f885c34ccc2679daa2c0007",
}
SEEDS = (314159, 271828)


def prepare(work):
    prepared = work / "prepared"
    prepared.mkdir(exist_ok=True, parents=True)
    sites = {}
    for code, digest in HASHES.items():
        path = work / f"{code}.pdb"
        if not path.exists():
            with urlopen(
                f"https://files.rcsb.org/download/{code}.pdb", timeout=60
            ) as response:
                path.write_bytes(response.read())
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Crystal source changed: {code}")
        lines = path.read_text().splitlines()
        if code == "1V9E":
            # The two independent polymer copies differ in the first coordination sphere.
            # B has three Zn--N contacts; A has two N contacts and a 2.27 Å CD2 contact.
            # Remap B's labels to the human-numbered A templates understood by receptor().
            lines = [
                line[:21] + "A" + f"{int(line[22:26]) + 1:4d}" + line[26:]
                for line in lines
                if line.startswith(("ATOM  ", "HETATM")) and line[21] == "B"
            ]
        sites[code] = receptor(prepared, code, lines, 20)
    ligands = list(csv.DictReader((HERE / "fluorinated.csv").open()))
    for ligand in ligands:
        name = ligand["ligand"]
        mol = Chem.MolFromSmiles(ligand["neutral_smiles"])
        ((donor,),) = mol.GetSubstructMatches(Chem.MolFromSmarts("[N;H2;$(NS(=O)=O)]"))
        mol.GetAtomWithIdx(donor).SetFormalCharge(-1)
        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol)
        opts = AllChem.ETKDGv3()
        opts.randomSeed = 2026
        if AllChem.EmbedMolecule(mol, opts) or AllChem.MMFFOptimizeMolecule(
            mol, maxIters=1000
        ):
            raise ValueError(name)
        (setup,) = MoleculePreparation().prepare(mol)
        text, ok, reason = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
        if not ok:
            raise ValueError(reason)
        serial = {}
        for line in text.splitlines():
            if line.startswith("REMARK INDEX MAP "):
                ids = list(map(int, line.split()[3:]))
                serial.update(zip(ids[::2], ids[1::2], strict=True))
        (prepared / f"{name}.pdbqt").write_text(text)
        for code, site in sites.items():
            point = " ".join(f"{v:.8f}" for v in site["target"])
            (prepared / f"{code}_{name}.anchor").write_text(
                f"{serial[donor + 1]} {point}\n"
            )
    (prepared / "sites.json").write_text(json.dumps(sites, indent=2) + "\n")
    return [r["ligand"] for r in ligands]


def dock(entry, work, binary):
    code, ligand, seed, arm = entry
    prepared = work / "prepared"
    prefix = work / "runs" / f"{code}_{ligand}_{seed}_{arm}"
    cmd = [
        str(binary),
        "--receptor",
        str(prepared / f"{code}_receptor.pdbqt"),
        "--config",
        str(prepared / f"{code}.box"),
        "--ligand",
        str(prepared / f"{ligand}.pdbqt"),
        "--seed",
        str(seed),
        "--cpu",
        "1",
        "--pose_selection",
        "late",
        "--exhaustiveness",
        "4",
        "--max_evals",
        "50000",
        "--energy_range",
        "100",
        "--num_modes",
        "9",
        "--search_counters",
        "--out",
        str(prefix.with_suffix(".pdbqt")),
    ]
    if arm == "point":
        cmd += ["--search_anchor", str(prepared / f"{code}_{ligand}.anchor")]
    with prefix.with_suffix(".log").open("w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", type=int, default=4)
    args = p.parse_args()
    work = args.work.resolve()
    names = prepare(work)
    (work / "runs").mkdir(exist_ok=True)
    with ThreadPoolExecutor(args.jobs) as pool:
        list(
            pool.map(
                lambda run: dock(run, work, args.binary.resolve()),
                itertools.product(HASHES, names, SEEDS, ("native", "point")),
            )
        )
