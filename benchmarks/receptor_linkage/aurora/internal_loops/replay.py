"""Fixed-coordinate checks and budget-matched blind cross-docking on 2X6E."""

import argparse
import hashlib
import json
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import geometry as geom
from dock import ARCHIVE
from pose_probe import mapped_rmsd

PROTEINS = {
    "dry": "prepared/2X6E.pdbqt",
    "1OL5_restrained": "template1OL5/trans_model/restrained/receptor.pdbqt",
    "1OL7_restrained": "template1OL7/trans_model/restrained/receptor.pdbqt",
    "1OL7_released": "template1OL7/trans_model/released/receptor.pdbqt",
}
ARMS = (*PROTEINS, "official", "double")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def invoke(args, log):
    start = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    with log.open("w") as f:
        subprocess.run(
            list(map(str, args)), stdout=f, stderr=subprocess.STDOUT, check=True
        )
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "wall_s": time.perf_counter() - start,
        "cpu_s": after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
    }


def parameters(work, drug):
    ligand = work / "prepared" / f"{drug}.pdbqt"
    bonds = (
        work
        / "raw"
        / ("at9283_chemical.bonds" if drug == "at9283" else "danusertib_native.bonds")
    )
    mol = Chem.MolFromMolFile(str(ligand.with_suffix(".sdf")), removeHs=True)
    (ref, comp, _) = geom.REFERENCES[drug][0]
    crystal = (work / (ref + ".pdb")).read_text().splitlines()
    ref_xyz, _ = geom.reference_mappings(
        mol, geom.ligand_atoms(crystal, comp, "A"), geom.ccd_data(work, comp), crystal
    )
    R, t, _ = geom.transform(
        geom.receiver_frames(ARCHIVE)["2X6E"], geom.ca_atoms(crystal, "A")
    )
    serials, _ = geom.serial_map(ligand, mol)
    return ligand, bonds, serials, ref_xyz @ R + t


def fixed(work, pro, out, arms):
    raw = out / "fixed"
    raw.mkdir(exist_ok=True)
    rows = []
    for drug in geom.REFERENCES:
        ligand, bonds, serials, refs = parameters(work, drug)
        for arm in arms:
            receptor = work / PROTEINS[arm]
            added = [
                np.array(geom.xyz(l))
                for l in receptor.read_text().splitlines()
                if l.startswith(("ATOM", "HETATM"))
                and int(l[22:26]) in (*range(281, 292), 303, 304, 305)
                and l.split()[-1] not in ("HD", "H")
            ]
            witness = None
            for i, pose in enumerate(("witness", *geom.SEEDS), start=1):
                p = work / "coordinates" / f"2X6E_{drug}_2X6E_pocket_{i}.pdbqt"
                xyz = geom.mode_coords(p, serials)[0]
                log = raw / f"{drug}_{arm}_{pose}.log"
                cost = invoke(
                    [
                        pro,
                        "--receptor",
                        receptor,
                        "--ligand",
                        ligand,
                        "--ligand_bonds",
                        bonds,
                        "--cpu",
                        1,
                        *[
                            arg
                            for dim, x in zip(
                                "xyz",
                                json.loads((work / "prepared/2X6E.json").read_text())[
                                    "center"
                                ],
                                strict=True,
                            )
                            for arg in (f"--center_{dim}", x, f"--size_{dim}", 26)
                        ],
                        "--pose_coordinates",
                        p,
                        "--unbound_energy",
                        0,
                        "--score_only",
                    ],
                    log,
                )
                score = float(
                    re.search(
                        r"Estimated Free Energy of Binding\s*:\s*([-+.\d]+)",
                        log.read_text(),
                    )[1]
                )
                if witness is None:
                    witness = score
                rows.append(
                    dict(
                        drug=drug,
                        arm=arm,
                        pose=pose,
                        score=score,
                        decoy_minus_witness=round(score - witness, 3),
                        rmsd_A=mapped_rmsd(xyz, refs),
                        nearest_added_A=cdist(xyz, added).min() if added else "",
                        input_sha256=digest(p),
                        **cost,
                    )
                )
    geom.write_csv(out, "fixed.csv", rows)


def blind(work, pro, official, out, arms):
    if (
        digest(official)
        != "f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644"
    ):
        raise ValueError("Expected the official v1.2.7 executable")
    raw = out / "blind"
    raw.mkdir(exist_ok=True)
    cells, modes = [], []
    center = json.loads((work / "prepared/2X6E.json").read_text())["center"]
    box = [
        arg
        for dim, x in zip("xyz", center, strict=True)
        for arg in (f"--center_{dim}", x, f"--size_{dim}", 26)
    ]
    for drug in geom.REFERENCES:
        ligand, _, serials, refs = parameters(work, drug)
        for i, seed in enumerate(geom.SEEDS):
            for arm in arms[i % len(arms) :] + arms[: i % len(arms)]:
                receptor = work / PROTEINS.get(arm, PROTEINS["dry"])
                binary = official if arm == "official" else pro
                output = raw / f"{arm}_{drug}_{seed}.pdbqt"
                cmd = [
                    binary,
                    "--receptor",
                    receptor,
                    "--ligand",
                    ligand,
                    *box,
                    "--seed",
                    seed,
                    "--exhaustiveness",
                    2,
                    "--cpu",
                    1,
                    "--num_modes",
                    9,
                    "--energy_range",
                    100,
                    "--max_evals",
                    50000 if arm == "double" else 25000,
                    "--out",
                    output,
                ]
                if arm != "official":
                    cmd += ["--pose_selection", "late", "--search_counters"]
                log = output.with_suffix(".log")
                cost = invoke(cmd, log)
                scores = [
                    float(l.split()[3])
                    for l in output.read_text().splitlines()
                    if l.startswith("REMARK VINA RESULT:")
                ]
                err = [
                    float(mapped_rmsd(p, refs))
                    for p in geom.mode_coords(output, serials)
                ]
                key = {"drug": drug, "seed": seed, "arm": arm}
                for rank, (r, s) in enumerate(zip(err, scores, strict=True), start=1):
                    modes.append(dict(**key, rank=rank, score=s, rmsd_A=r))
                counts = re.findall(
                    r"\b(?:hunt|second|explicit)_evals=(\d+)", log.read_text()
                )
                cells.append(
                    dict(
                        **key,
                        n=len(err),
                        top1=err[0],
                        best5=min(err[:5]),
                        oracle=min(err),
                        regret=err[0] - min(err),
                        evals=sum(map(int, counts)) if counts else "",
                        binary_sha256=digest(binary),
                        receptor_sha256=digest(receptor),
                        pool_sha256=digest(output),
                        **cost,
                    )
                )
    geom.write_csv(out, "modes.csv", modes)
    geom.write_csv(out, "cells.csv", cells)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=("blind", "fixed"))
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("pro", type=Path)
    p.add_argument("--official", type=Path)
    p.add_argument("--arms", choices=ARMS, nargs="+", default=ARMS)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    if a.stage == "blind":
        blind(
            a.work.resolve(),
            a.pro.resolve(),
            a.official.resolve(),
            a.out.resolve(),
            a.arms,
        )
    else:
        fixed(
            a.work.resolve(),
            a.pro.resolve(),
            a.out.resolve(),
            [k for k in a.arms if k in PROTEINS],
        )
