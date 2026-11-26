"""Replay the human PNP Table 1 rate-based inhibition panel."""

import argparse
import csv
import itertools
import json
import statistics
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from openmm import Platform
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem

sys.path.insert(
    0, str(Path(__file__).resolve().parents[1] / "kinase_state_affinity" / "source")
)
from prepare import ligand

SEEDS = (271828, 314159, 161803, 141421)
FLAGS = {
    "native": [],
    "balanced": ["--balanced_torsion_search"],
    "prepared": ["--prepared_torsion_starts"],
}
HERE = Path(__file__).resolve().parent


def prepare(work, table):
    p = work / "prepared"
    p.mkdir(exist_ok=True)
    rows = (work / "raw/1RR6.pdb").read_text().splitlines()
    chains = [[], [], []]
    ions = []
    bound = []
    for line in rows:
        if line.startswith("MODEL"):
            model = int(line.split()[1]) - 1
        if (
            line[:6] not in ("HETATM", "ATOM  ")
            or line[16] not in " A"
            or line[76:78].strip() in ("H", "D")
        ):
            continue
        line = line[:16] + " " + line[17:21] + chr(65 + model) + line[22:]
        if line.startswith("ATOM"):
            chains[model].append(line)
        elif line[17:20] == "PO4":
            ions.append(line)
        elif model == 0 and line[17:20] == "IMH":
            bound.append([float(line[k : k + 8]) for k in (30, 38, 46)])
    assert len(bound) == 19 and all(chains) and ions
    (p / "box.txt").write_text(
        "".join(
            f"center_{q} = {v:.5f}\nsize_{q} = 22\n"
            for q, v in zip("xyz", np.mean(bound, axis=0))
        )
    )
    inp = p / "protein.pdb"
    inp.write_text("\nTER\n".join("\n".join(chain) for chain in chains) + "\nEND\n")
    f = PDBFixer(filename=str(inp), platform=Platform.getPlatformByName("CPU"))
    f.findMissingResidues()
    f.missingResidues = {}
    f.findMissingAtoms()
    missing = {
        f"{r.chain.id}:{r.name}{r.id}": [a.name for a in atoms]
        for r, atoms in f.missingAtoms.items()
    }
    (p / "missing.json").write_text(json.dumps(missing, indent=2) + "\n")
    f.addMissingAtoms(seed=2026)
    with inp.open("w") as out:
        PDBFile.writeFile(f.topology, f.positions, out, keepIds=True)
    with (p / "receptor.log").open("w") as log:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "meeko.cli.mk_prepare_receptor",
                "--read_pdb",
                inp,
                "-p",
                p / "protein.pdbqt",
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    protein = (p / "protein.pdbqt").read_text()
    index = max(
        int(l[6:11]) for l in protein.splitlines() if l.startswith(("ATOM", "HETATM"))
    )
    with (p / "receptor.pdbqt").open("w") as out:
        out.write(protein)
        for i, l in enumerate(ions, index + 1):
            elem = l[76:78].strip()
            # Acceptor-only inorganic phosphate approximation. Native Vina does not use partial charges.
            q, typ = (1.0, "P") if elem == "P" else (-0.75, "OA")
            out.write(l[:6] + f"{i:5d}" + l[11:66] + f"    {q:6.3f} {typ:2s}\n")
    for r in table:
        smi = json.loads((work / "raw" / (r["ccd"] + ".json")).read_text())[
            "rcsb_chem_comp_descriptor"
        ]["SMILES_stereo"]
        if r["ligand"].endswith("ImmG"):
            assert smi.count("NC=N2") == 1
            smi = smi.replace("NC=N2", "NC(N)=N2")
        mol = Chem.MolFromSmiles(smi)
        if Chem.GetFormalCharge(mol) == 0:
            [a] = [
                a
                for a in mol.GetAtoms()
                if a.GetAtomicNum() == 7
                and not a.GetIsAromatic()
                and all(
                    n.GetAtomicNum() == 6 and not n.GetIsAromatic()
                    for n in a.GetNeighbors()
                )
            ]
            a.SetNumExplicitHs(a.GetTotalNumHs() + 1)
            a.SetFormalCharge(1)
            a.SetNoImplicit(True)
            Chem.SanitizeMol(mol)
        assert Chem.GetFormalCharge(mol) == 1
        smi = Chem.MolToSmiles(mol)
        (p / (r["ligand"] + ".smi")).write_text(smi + "\n")
        ligand(work, {"cid": r["ligand"], "smiles": smi})


def run(work, binary, table, jobs):
    dest = work / "pools"
    dest.mkdir(exist_ok=True)

    def one(case):
        name, seed, arm = case
        base = dest / f"{name}_{seed}_{arm}"
        with base.with_suffix(".log").open("w") as log:
            subprocess.run(
                [
                    binary,
                    "--receptor",
                    work / "prepared/receptor.pdbqt",
                    "--ligand",
                    work / "prepared" / (name + ".pdbqt"),
                    "--config",
                    work / "prepared/box.txt",
                    "--pose_selection",
                    "late",
                    "--seed",
                    str(seed),
                    "--cpu",
                    "1",
                    "--exhaustiveness",
                    "2",
                    "--max_evals",
                    "25000",
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--out",
                    base.with_suffix(".pdbqt"),
                    *FLAGS[arm],
                ],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
            )

    with ThreadPoolExecutor(max_workers=jobs) as e:
        list(e.map(one, itertools.product([r["ligand"] for r in table], SEEDS, FLAGS)))


def summarize(work, table):
    results = []
    pairs = []
    for endpoint in ("Ki25", "KiSlow25", "Ki37", "KiSlow37"):
        for scope in ("all", "exclude_cocrystal_ImmH"):
            rows = [r for r in table if scope == "all" or r["ligand"] != "ImmH"]
            for tag, seeds in (
                ("all", SEEDS),
                ("first", SEEDS[:2]),
                ("second", SEEDS[2:]),
            ):
                scores = {}
                for r, arm in itertools.product(rows, FLAGS):
                    E = []
                    for s in seeds:
                        p = work / "pools" / f"{r['ligand']}_{s}_{arm}.pdbqt"
                        E.append(
                            float(
                                next(
                                    l
                                    for l in p.read_text().splitlines()
                                    if l.startswith("REMARK VINA RESULT")
                                ).split()[3]
                            )
                        )
                    scores[r["ligand"], arm] = statistics.median(E)
                edges = [
                    (a, b)
                    for a, b in itertools.permutations(rows, 2)
                    if float(b[endpoint]) >= 3 * float(a[endpoint])
                ]
                for arm in FLAGS:
                    cuts = []
                    for a, b in edges:
                        strong, weak = a["ligand"], b["ligand"]
                        gap = scores[weak, arm] - scores[strong, arm]
                        correct = 1 if gap > 1e-9 else 0 if gap < -1e-9 else 0.5
                        cuts.append(correct)
                        if tag == "all" and scope == "all":
                            pairs.append(
                                [
                                    endpoint,
                                    arm,
                                    strong,
                                    weak,
                                    a[endpoint],
                                    b[endpoint],
                                    gap,
                                    correct,
                                ]
                            )
                    results.append([endpoint, scope, tag, arm, len(edges), sum(cuts)])
    for name, keys, data in [
        (
            "summary",
            ["endpoint", "scope", "seed_group", "arm", "n", "correct"],
            results,
        ),
        (
            "pairs",
            [
                "endpoint",
                "arm",
                "strong",
                "weak",
                "stronger_pM",
                "weaker_pM",
                "predicted_gap",
                "correct",
            ],
            pairs,
        ),
    ]:
        with (work / (name + ".csv")).open("w") as out:
            w = csv.writer(out)
            w.writerow(keys)
            w.writerows(data)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", type=int, default=3)
    args = p.parse_args()
    table = list(csv.DictReader((HERE / "table1.csv").open()))
    prepare(args.work, table)
    run(args.work, args.binary.resolve(), table, args.jobs)
    summarize(args.work, table)
