"""Independent cocrystal/ITC panel from Baum et al. JMB 2009, DOI 10.1016/j.jmb.2009.04.051."""

import argparse
import csv
import itertools
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

# Compound, P1 ring attached to benzyl CH2, original ITC DeltaG (kJ/mol), cocrystal if available.
# deltaG recorded from the measured Kd in paper Table 1. All assay temperatures are 298.15 K.
SERIES = {
    "4": ("c1ccc(C(=[NH2+])N)cc1", -46.1, "2ZDA", "32U"),
    "6a": ("c1cc(Cl)c(Cl)cc1", -38.4, None, None),
    "6e": ("c1cc(Cl)c(F)cc1", -37.3, None, None),
    "6b": ("c1cc(Cl)c(C)cc1", -37.2, None, None),
    "1c": ("c1cccc(Br)c1", -35.8, None, None),
    "1b": ("c1cccc(Cl)c1", -35.4, "2ZC9", "22U"),
    "3a": ("c1cccc(C)c1", -34.8, "2ZF0", "51U"),
    "1d": ("c1cccc(I)c1", -34.5, None, None),
    "7a": ("c1cc(C)c(C)cc1", -34.4, None, None),
    "1a": ("c1cccc(F)c1", -31.3, "2ZDV", "37U"),
    "3b": ("c1cccc(CC)c1", -32.9, None, None),
    "5": ("c1ccccc1", -31.7, "2ZFF", "53U"),
}
CORE = "[NH3+][C@H](Cc1ccccc1)C(=O)N1CCC[C@H]1C(=O)NC"
REC = ("3D49", "2ZFF", "2ZC9")
SEEDS = (83519, 194027)
POCKET = (
    (57, " "),
    (60, "A"),
    (60, "D"),
    (99, " "),
    (174, " "),
    (189, " "),
    (195, " "),
    (215, " "),
    (216, " "),
    (227, " "),
    (228, " "),
)


def xyz(r):
    return np.array([float(r[i : i + 8]) for i in (30, 38, 46)])


def lines(p):
    return [
        r
        for r in Path("raw", p + ".pdb").read_text().splitlines()
        if r.startswith(("ATOM  ", "HETATM"))
        and r[16] in " A"
        and r[76:78].strip() not in ("H", "D")
    ]


def ca(p):
    coordinates = {
        (int(r[22:26]), r[26]): xyz(r)
        for r in lines(p)
        if r[21] == "H" and r[12:16] == " CA " and r.startswith("ATOM")
    }
    return np.array([coordinates[k] for k in POCKET])


def transform(p, rec):
    a, b = ca(p), ca(rec)
    u, _s, vh = np.linalg.svd((a - a.mean(0)).T @ (b - b.mean(0)))
    u[:, -1] *= np.linalg.det(u @ vh)
    rot = u @ vh
    return rot, b.mean(0) - a.mean(0) @ rot


def confs(ident, ring):
    mol = Chem.AddHs(Chem.MolFromSmiles(CORE + ring))
    assert [x[1] for x in Chem.FindMolChiralCenters(mol)] == ["R", "S"]
    opts = AllChem.ETKDGv3()
    opts.randomSeed = 408151
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=8, params=opts)
    vals = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=1200, numThreads=1)
    fit = [(e, cid) for cid, (status, e) in zip(ids, vals) if not status]
    if not fit:
        raise RuntimeError(ident)
    best = Chem.Mol(mol, confId=min(fit)[1])
    (st,) = MoleculePreparation().prepare(best)
    text, ok, why = PDBQTWriterLegacy.write_string(st, add_index_map=True)
    if not ok:
        raise ValueError(why)
    Path("prepared", ident + ".pdbqt").write_text(text)


def prepare():
    for path in ("prepared", "results"):
        Path(path).mkdir(exist_ok=True)
    for ident, (ring, *_) in SERIES.items():
        confs(ident, ring)
    for rec in REC:
        rec_lines = [
            r[:16] + " " + r[17:]
            for r in lines(rec)
            if r.startswith("ATOM  ") and r[21] in "HL"
        ]
        inp = Path("prepared", rec + ".pdb")
        inp.write_text("\n".join(rec_lines) + "\nEND\n")
        with Path("prepared", rec + ".log").open("w") as out:
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "meeko.cli.mk_prepare_receptor",
                    "--read_pdb",
                    str(inp),
                    "-a",
                    "-p",
                    str(inp.with_suffix(".pdbqt")),
                ],
                stdout=out,
                stderr=subprocess.STDOUT,
                check=True,
            )
        # Box from common binding pocket protein landmarks only, no query ligand coordinates.
        center = ca(rec).mean(0)
        Path("prepared", rec + ".box").write_text(
            "".join(
                f"center_{ax} = {x:.5f}\nsize_{ax} = 25\n"
                for ax, x in zip("xyz", center)
            )
        )
    with Path("series.csv").open("w") as h:
        w = csv.writer(h)
        w.writerow(
            [
                "id",
                "smiles",
                "ITC_deltaG_kJ_per_mol",
                "Kd_uM_derived_at_298.15K",
                "PDB",
                "CCD",
            ]
        )
        for id, (ring, dg, p, c) in SERIES.items():
            w.writerow(
                [
                    id,
                    CORE + ring,
                    dg,
                    1e6 * np.exp(dg / (0.00831446261815324 * 298.15)),
                    p,
                    c,
                ]
            )


def dock(cell, binary):
    rec, ident, seed = cell
    basename = f"{rec}_{ident}_{seed}"
    outf = Path("results", basename + ".pdbqt")
    with Path("results", basename + ".log").open("w") as log:
        subprocess.run(
            [
                binary,
                "--config",
                f"prepared/{rec}.box",
                "--receptor",
                f"prepared/{rec}.pdbqt",
                "--ligand",
                f"prepared/{ident}.pdbqt",
                "--cpu",
                "1",
                "--pose_selection",
                "late",
                "--exhaustiveness",
                "4",
                "--max_evals",
                "50000",
                "--num_modes",
                "12",
                "--refine_candidates",
                "16",
                "--energy_range",
                "8",
                "--seed",
                str(seed),
                "--out",
                str(outf),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("binary")
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--prepared", action="store_true")
    args = ap.parse_args()
    if not args.prepared:
        prepare()
    Path("results").mkdir(exist_ok=True)
    # Run the five cocrystal compounds first, followed by the remaining ITC compounds.
    ids = sorted(SERIES, key=lambda i: SERIES[i][2] is None)
    cells = [(rec, i, seed) for i, seed, rec in itertools.product(ids, SEEDS, REC)]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(lambda cell: dock(cell, args.binary), cells))
