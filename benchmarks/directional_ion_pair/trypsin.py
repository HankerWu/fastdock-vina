"""Independent bovine trypsin crystal and same-buffer calorimetric transfers."""

import argparse
import itertools
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

from experiment import ARMS, SEEDS, aspsite, prepare, run, write

# Talhout & Engberts, Eur J Biochem 2001, Table 1 and parental benzamidine at 25 C.
# para-CONH2 is the formula in Tables 1 and 2. The adjacent prose calls it amidino.
SERIES = {
    "H": ("NC(=[NH2+])c1ccccc1", 4.50),
    "NH2": ("NC(=[NH2+])c1ccc(N)cc1", 12.7),
    "CH3": ("NC(=[NH2+])c1ccc(C)cc1", 6.94),
    "OCH3": ("NC(=[NH2+])c1ccc(OC)cc1", 2.65),
    "CONH2": ("NC(=[NH2+])c1ccc(C(N)=O)cc1", 1.52),
}
RECEPTORS = ("3PTB", "1BTY")


def prep(out, raw):
    raw.mkdir(exist_ok=True, parents=True)
    for rec in RECEPTORS:
        path = raw / f"{rec}.pdb"
        if not path.exists():
            path.write_bytes(
                urllib.request.urlopen(
                    f"https://files.rcsb.org/download/{rec}.pdb", timeout=40
                ).read()
            )
        lines = path.read_text().splitlines()
        heavy = [
            l[:16] + " " + l[17:]
            for l in lines
            if l.startswith("ATOM  ")
            and l[16] in " A"
            and l[21] == "A"
            and l[76:78].strip() not in ("H", "D")
        ]
        (out / f"{rec}.protein.pdb").write_text("\n".join(heavy) + "\nEND\n")
        with (out / f"{rec}.prep.log").open("w") as log:
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "meeko.cli.mk_prepare_receptor",
                    "--read_pdb",
                    out / f"{rec}.protein.pdb",
                    "-p",
                    out / f"{rec}.pdbqt",
                ],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        cg, o1, o2 = aspsite(path, "A")
        normal = (o1 + o2) / 2 - cg
        center = cg + 6 * normal / np.linalg.norm(normal)
        (out / f"{rec}.box").write_text(
            "".join(
                f"center_{a} = {v:.4f}\nsize_{a} = 26\n" for a, v in zip("xyz", center)
            )
        )
    for label, (smiles, ka) in SERIES.items():
        m = Chem.AddHs(Chem.MolFromSmiles(smiles))
        par = AllChem.ETKDGv3()
        par.randomSeed = 2026
        if AllChem.EmbedMolecule(m, par) or AllChem.MMFFOptimizeMolecule(
            m, maxIters=1000
        ):
            raise ValueError(label)
        (setup,) = MoleculePreparation().prepare(m)
        text, ok, why = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
        if not ok:
            raise ValueError(why)
        (out / f"{label}.pdbqt").write_text(text)
    for rec in RECEPTORS:
        prepare(out, rec, raw / f"{rec}.pdb", "A", [out / f"{k}.pdbqt" for k in SERIES])
    write(
        out / "series.csv",
        [
            dict(ligand=k, smiles=s, ka_M_inverse=v * 1e4, kd_uM=100 / v)
            for k, (s, v) in SERIES.items()
        ],
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("raw", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    p.add_argument("--arms", nargs="+", choices=ARMS, default=ARMS)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    prep(a.out, a.raw)
    inputs = {
        (r, l): (a.out / f"{l}.pdbqt", a.out / f"{r}.pdbqt", a.out / f"{r}.box")
        for r in RECEPTORS
        for l in SERIES
    }
    run(
        a.binary.resolve(),
        a.out,
        inputs,
        itertools.product(RECEPTORS, SERIES, a.seeds, a.arms),
        a.jobs,
    )
