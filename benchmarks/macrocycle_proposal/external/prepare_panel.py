"""Prepare prospective cognate macrocycles from non-crystal CCD-derived coordinates."""

import argparse
import json
import os
import sys
import time
from pathlib import Path

CASES = {
    "4CLI": ("ALK", "5P8", "A"),
    "8BI2": ("SYK", "QR7", "A"),
    "8K78": ("MET", "IYC", "A"),
    "1YND": ("PPIA", "SFA", "A"),
    "7MBO": ("F11", "YXG", "A"),
    "5QQP": ("F11", "NR7", "A"),
    "5T9W": ("PPIA", "78E", "A"),
    "1NMK": ("PPIA", "SFM", "A"),
    "8K79": ("SRC", "IYC", "A"),
    "6FS0": ("MCL1", "E4W", "A"),
}


def run(repo: Path, work: Path, accession: str):
    import gemmi
    from rdkit import Chem

    sys.path.insert(0, str(repo / "benchmarks/pair_gradient/pose"))
    import prepare

    prepare.CASES.update(CASES)
    started = time.perf_counter()
    directory = work / accession
    directory.mkdir(exist_ok=True, parents=True)
    record = {
        "pdb": accession,
        "target": CASES[accession][0],
        "ccd": CASES[accession][1],
        "chain": "A",
    }
    try:
        prepare.prepare(work, accession)
        prepared = json.loads((directory / "preparation.json").read_text())
        mol = Chem.MolFromSmiles(prepared["canonical_smiles"])
        ring_sizes = [len(ring) for ring in mol.GetRingInfo().AtomRings()]
        ligand = (directory / "ligand.pdbqt").read_text().splitlines()
        types = [l.split()[-1] for l in ligand if l.startswith(("ATOM ", "HETATM"))]
        status = prepared["receptor_status"]
        if status == "ready" and not any(t.startswith("CG") for t in types):
            status = "no opened macrocycle in default Meeko PDBQT"
        if status == "ready" and accession == "6FS0":
            status = "excluded: AZD5991 atropisomer not constrained by the CCD SMILES"
        pdb = gemmi.read_structure(str(work / "raw" / (accession + ".pdb")))
        record.update(
            status=status,
            resolution_A=pdb.resolution,
            ring_sizes=ring_sizes,
            atoms=mol.GetNumHeavyAtoms(),
            formal_charge=Chem.GetFormalCharge(mol),
            glue_labels=[t for t in types if t.startswith(("CG", "G"))],
            maps=len(prepared["mapped"]),
            pdb_sha256=prepared["pdb_sha256"],
            ccd_sha256=prepared["ccd_sha256"],
            nearest_omitted=prepared["closest_omitted_heavy_atom"],
            nearest_repaired=prepared.get("closest_repaired_ligand_distance"),
            torsdof=int(next(l for l in ligand if l.startswith("TORSDOF")).split()[1]),
        )
    except (ValueError, KeyError, IndexError, RuntimeError) as exc:
        record["status"] = str(exc)
    record["preparation_wall_s"] = time.perf_counter() - started
    (directory / "screening.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("repo", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("ids", nargs="*")
    args = p.parse_args()
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    for accession in args.ids or CASES:
        run(args.repo.resolve(), args.work.resolve(), accession)
