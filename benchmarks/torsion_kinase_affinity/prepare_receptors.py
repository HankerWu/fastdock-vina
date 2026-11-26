"""Complete distal missing side chains conservatively, and preserve observed pocket heavy atoms."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

os.environ.setdefault("OPENMM_CPU_THREADS", "1")
import numpy as np
from openmm import Platform
from openmm.app import PDBFile
from pdbfixer import PDBFixer


def xyz(line):
    return np.array([float(line[k : k + 8]) for k in (30, 38, 46)])


def prepare(work, pdb, component):
    raw = (work / "raw" / f"{pdb}.pdb").read_text().splitlines()
    out = work / "prepared"
    out.mkdir(exist_ok=True)
    lig = np.asarray(
        [
            xyz(l)
            for l in raw
            if l.startswith("HETATM")
            and l[21] == "A"
            and l[17:20] == component
            and l[16] in " A"
            and l[76:78].strip() not in ("H", "D")
        ]
    )
    if len(lig) < 9:
        raise ValueError(component)
    near = lambda points: float(
        np.linalg.norm(np.asarray(points)[:, None, :] - lig[None, :, :], axis=2).min()
    )
    records = [
        l[:16] + " " + l[17:]
        for l in raw
        if l.startswith("ATOM  ")
        and l[21] == "A"
        and l[16] in " A"
        and l[76:78].strip() not in ("H", "D")
    ]
    file = out / (pdb + "-protein.pdb")
    file.write_text("\n".join(records) + "\nEND\n")
    f = PDBFixer(filename=str(file), platform=Platform.getPlatformByName("CPU"))
    f.findMissingResidues()
    f.missingResidues = {}
    f.findMissingAtoms()
    truncations = {}
    for r, atoms in f.missingAtoms.items():
        obs = [l for l in records if l[22:26].strip() == r.id]
        dist = near([xyz(l) for l in obs])
        # Retain observed backbone and CB of a distal, incomplete residue.
        if dist > 15 and {"N", "CA", "C", "O", "CB"} <= {l[12:16].strip() for l in obs}:
            truncations[r.id] = dict(
                original=r.name, distance=dist, unresolved=[a.name for a in atoms]
            )
    clean = []
    for l in records:
        if l[22:26].strip() in truncations:
            if l[12:16].strip() not in {"N", "CA", "C", "O", "CB", "OXT"}:
                continue
            l = l[:17] + "ALA" + l[20:]
        clean.append(l)
    file.write_text("\n".join(clean) + "\nEND\n")
    f = PDBFixer(filename=str(file), platform=Platform.getPlatformByName("CPU"))
    f.findMissingResidues()
    gaps = {str(k): v for k, v in f.missingResidues.items()}
    f.missingResidues = {}
    f.findMissingAtoms()
    missing = {r.id: [a.name for a in atoms] for r, atoms in f.missingAtoms.items()}
    f.addMissingAtoms(seed=43103)
    with file.open("w") as fd:
        PDBFile.writeFile(f.topology, f.positions, fd, keepIds=True)
    repair = [
        xyz(l)
        for l in file.read_text().splitlines()
        if l.startswith("ATOM  ")
        and l[12:16].strip() in missing.get(l[22:26].strip(), [])
    ]
    close = near(repair) if repair else None
    (out / (pdb + ".json")).write_text(
        json.dumps(
            dict(
                truncations=truncations,
                missing=missing,
                closest_repaired=close,
                reported_gaps=gaps,
                center=lig.mean(axis=0).tolist(),
            ),
            indent=2,
        )
        + "\n"
    )
    with (out / (pdb + ".prep.log")).open("w") as fd:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "meeko.cli.mk_prepare_receptor",
                "--read_pdb",
                file,
                "-p",
                out / (pdb + ".pdbqt"),
            ],
            stdout=fd,
            stderr=subprocess.STDOUT,
            check=True,
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    args = p.parse_args()
    for pdb, component in [("2X6E", "YM4"), ("2J4Z", "626")]:
        prepare(args.work, pdb, component)
