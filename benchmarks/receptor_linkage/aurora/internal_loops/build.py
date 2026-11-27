"""Restore only the two unresolved internal sequences to the frozen 2X6E scaffold."""

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENMM_CPU_THREADS", "1")
import numpy as np
from openmm import Platform, app
from pdbfixer import PDBFixer

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "torsion_kinase_affinity"))
from prepare_receptors import deposited_gaps, xyz


def lines(path):
    return [
        r for r in path.read_text().splitlines() if r.startswith(("ATOM  ", "HETATM"))
    ]


def atom_key(row):
    return (row[21], int(row[22:26]), row[12:16].strip())


def model(work, seed, repaired=False):
    start = time.perf_counter()
    prep = work / "prepared"
    out = work / f"loop{seed}"
    out.mkdir(exist_ok=True)
    archive_pdb = prep / "2X6E-protein.pdb"
    pdb_path = out / ("repaired.pdb" if repaired else "completed.pdb")
    if not repaired:
        fixer = PDBFixer(
            filename=str(archive_pdb), platform=Platform.getPlatformByName("CPU")
        )
        residues = list(fixer.topology.residues())
        # Map by deposited flanking residue IDs to the archived mutated (distal Ala) scaffold.
        raw_path = work / "raw" / "2X6E.pdb"
        raw = app.PDBFile(str(raw_path))
        source_chains = list(raw.topology.chains())
        source_index = {}
        for (chain, pos), names in deposited_gaps(raw_path, "A").items():
            original = list(source_chains[chain].residues())
            if 0 < pos < len(original):
                source_index[int(original[pos].id)] = names
        insertion = {
            (0, i): source_index[int(r.id)]
            for i, r in enumerate(residues)
            if int(r.id) in source_index
        }
        if [len(v) for v in insertion.values()] != [11, 3]:
            raise ValueError(insertion)
        fixer.missingResidues = insertion
        fixer.findMissingAtoms()
        if any(fixer.missingAtoms.values()):
            raise ValueError("The archived scaffold unexpectedly lacks heavy atoms")
        fixer.addMissingAtoms(seed=seed)
        with pdb_path.open("w") as f:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    before = {atom_key(r): xyz(r) for r in lines(archive_pdb)}
    after = {atom_key(r): xyz(r) for r in lines(pdb_path)}
    displacement = max(np.linalg.norm(after[k] - before[k]) for k in before)
    if displacement > 0.002:
        raise ValueError(("original atom moved", displacement))
    new = {k for k in after if k not in before}
    new_res = {i for _, i, _ in new}
    if new_res != set(range(281, 292)) | {303, 304, 305}:
        raise ValueError(new_res)
    peptide = [
        float(np.linalg.norm(after["A", i, "C"] - after["A", i + 1, "N"]))
        for i in sorted(new_res | {280, 302})
    ]
    chirality = {
        i: float(
            np.linalg.det(
                [after["A", i, a] - after["A", i, "CA"] for a in ("N", "C", "CB")]
            )
        )
        for i in sorted(new_res)
        if ("A", i, "CB") in after
    }
    # Exclude covalently adjacent residues when looking for severe packing overlaps.
    near_nonlocal = min(
        float(np.linalg.norm(after[k] - p))
        for k in new
        for b, p in after.items()
        if abs(k[1] - b[1]) > 1
    )
    data = {
        "seed": seed,
        "postrepair": repaired,
        "added_heavy": len(new),
        "max_original_A": displacement,
        "rebuilt_peptide_A": peptide,
        "signed_N_C_CB_A3": chirality,
        "closest_nonadjacent_loop_contact_A": near_nonlocal,
        "elapsed_s": time.perf_counter() - start,
    }
    label = "repaired" if repaired else "raw"
    (out / f"{label}_geometry.json").write_text(json.dumps(data, indent=2) + "\n")
    return data


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("seed", type=int)
    p.add_argument("--repaired", action="store_true")
    args = p.parse_args()
    print(json.dumps(model(args.work.resolve(), args.seed, args.repaired)))
