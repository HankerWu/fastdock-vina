"""Prepare a receptor and an independently embedded ligand, then sample native poses."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import cdist

SEEDS = (314159, 271828, 161803)


def prepare(source, folder, case):
    folder.mkdir(exist_ok=True, parents=True)
    for old, new in [("protein.pdb", "protein.pdb"), ("ligand.sdf", "crystal.sdf")]:
        shutil.copyfile(source / case / (case + "_" + old), folder / new)
    ref = Chem.SDMolSupplier(str(folder / "crystal.sdf"), removeHs=True)[0]
    xyz = ref.GetConformer().GetPositions()
    # Molecular graphs and stereochemistry survive, the deposited coordinates do not.
    query = Chem.Mol(ref)
    query.RemoveAllConformers()
    query = Chem.AddHs(query)
    params = AllChem.ETKDGv3()
    params.randomSeed = 150473
    if AllChem.EmbedMolecule(query, params) < 0:
        raise ValueError("ligand embed failed")
    if AllChem.MMFFOptimizeMolecule(query, maxIters=2000) != 0:
        raise ValueError("MMFF optimization failed")
    (setup,) = MoleculePreparation().prepare(query)
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(error)
    (folder / "query.pdbqt").write_text(pdbqt)
    (folder / "box").write_text(
        "".join(
            f"center_{a} = {v:.5f}\nsize_{a} = 26\n" for a, v in zip("xyz", xyz.mean(0))
        )
    )
    protein = [
        line
        for line in (folder / "protein.pdb").read_text().splitlines()
        if line.startswith("ATOM  ")
        and line[16] in " A"
        and line[76:78].strip() not in ("H", "D")
    ]
    (folder / "polymer.pdb").write_text("\n".join(protein) + "\nEND\n")
    with (folder / "prepare.log").open("w") as fd:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "meeko.cli.mk_prepare_receptor",
                "--read_pdb",
                str(folder / "polymer.pdb"),
                "-a",
                "-p",
                str(folder / "receptor.pdbqt"),
            ],
            stdout=fd,
            stderr=subprocess.STDOUT,
            check=True,
        )
    # Reject omission of any observed heavy atom in a near-ligand polymer residue.
    spec = lambda l: (l[21], l[22:27], l[17:20], l[12:16].strip())
    keys = {
        spec(l)
        for l in (folder / "receptor.pdbqt").read_text().splitlines()
        if l.startswith(("ATOM  ", "HETATM"))
    }
    dropped = [l for l in protein if spec(l) not in keys]
    if dropped:
        distances = cdist(
            [[float(l[s : s + 8]) for s in (30, 38, 46)] for l in dropped], xyz
        ).min(axis=1)
        if min(distances) < 6:
            raise ValueError(
                f"near ligand missing atom ({min(distances):.2f} A): {dropped[np.argmin(distances)][12:27]}"
            )
    (folder / "ready").write_text("prepared\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument(
        "--redock", action="store_true", help="rerun completed pose pools"
    )
    args = parser.parse_args()
    cases = (Path(__file__).with_name("panel.txt")).read_text().split()

    def process(case):
        folder = args.work / case
        try:
            if not (folder / "ready").exists():
                prepare(args.source, folder, case)
            for seed in SEEDS:
                out = folder / f"poses-{seed}.pdbqt"
                if (
                    not args.redock
                    and out.exists()
                    and out.read_text().endswith("ENDMDL\n")
                ):
                    continue
                cmd = [
                    str(args.binary.resolve()),
                    "--receptor",
                    str(folder / "receptor.pdbqt"),
                    "--ligand",
                    str(folder / "query.pdbqt"),
                    "--config",
                    str(folder / "box"),
                    "--seed",
                    str(seed),
                    "--cpu",
                    "1",
                    "--exhaustiveness",
                    "3",
                    "--max_evals",
                    "60000",
                    "--num_modes",
                    "9",
                    "--pose_selection",
                    "late",
                    "--refine_candidates",
                    "12",
                    "--energy_range",
                    "100",
                    "--out",
                    str(out),
                ]
                with (folder / f"dock-{seed}.log").open("w") as fd:
                    subprocess.run(cmd, stdout=fd, stderr=subprocess.STDOUT, check=True)
            print(case, "finished", flush=True)
        except Exception as e:
            (folder / "excluded.txt").write_text(str(e) + "\n")
            print(case, str(e), flush=True)

    with ThreadPoolExecutor(args.jobs) as executor:
        list(executor.map(process, cases))


if __name__ == "__main__":
    main()
