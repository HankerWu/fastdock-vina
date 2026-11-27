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
from openmm.app import PDBFile, Topology
from pdbfixer import PDBFixer


def xyz(line):
    return np.array([float(line[k : k + 8]) for k in (30, 38, 46)])


def deposited_gaps(path, chain_id):
    """Read polymer sequence before ATOM filtering removes SEQRES records."""
    structure = PDBFixer(filename=str(path))
    structure.findMissingResidues()
    chains = list(structure.topology.chains())
    return {
        key: names for key, names in structure.missingResidues.items()
        if chains[key[0]].id == chain_id
    }


def write_observed_segments(path, topology, positions, gap_starts):
    """Preserve missing polymer segments in the coordinate export read by OpenMM."""
    segments = Topology()
    replacement = {}
    for old_chain in topology.chains():
        chain = segments.addChain(old_chain.id)
        for pos, residue in enumerate(old_chain.residues()):
            if pos in gap_starts and pos:
                chain = segments.addChain(old_chain.id)
            new = segments.addResidue(residue.name, chain, residue.id, residue.insertionCode)
            for a in residue.atoms():
                replacement[a] = segments.addAtom(a.name, a.element, new, a.id)
    for bond in topology.bonds():
        a, b = (replacement[atom] for atom in bond)
        if a.residue.chain != b.residue.chain and {a.name, b.name} == {"C", "N"}:
            continue
        segments.addBond(a, b, bond.type, bond.order)
    segments.setPeriodicBoxVectors(topology.getPeriodicBoxVectors())
    with path.open("w") as f:
        PDBFile.writeFile(segments, positions, f, keepIds=True)


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
    gaps = deposited_gaps(work / "raw" / f"{pdb}.pdb", "A")
    res = list(f.topology.residues())
    internal = {p: names for (_, p), names in gaps.items() if 0 < p < len(res)}
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
                reported_gaps={str(k): v for k, v in gaps.items()},
                unresolved_internal_segments=[
                    dict(after=res[p-1].id, before=res[p].id, sequence=names)
                    for p, names in internal.items()
                ],
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

    # The rigid PDBQT uses its established templates. The coordinate export must
    # also mark gaps so an all-atom reader cannot bond the two distant anchors.
    if internal:
        write_observed_segments(file, f.topology, f.positions, internal)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    args = p.parse_args()
    for pdb, component in [("2X6E", "YM4"), ("2J4Z", "626")]:
        prepare(args.work, pdb, component)
