"""Document every input preparation attempt, source, version and structural limitation."""

import argparse
import json
import sys
from importlib.metadata import version
from pathlib import Path
from urllib.request import urlopen

import gemmi

from external_run import sha256
from prepare_panel import CASES

BOUND = {"5T9W": 300, "7MBO": 302}


def main(repo, work):
    sys.path.insert(0, str(repo / "benchmarks/pair_gradient/pose"))
    from sources import crystal_contacts

    records = []
    for pdb, (target, ccd, chain) in CASES.items():
        d = work / pdb
        row = json.loads((d / "screening.json").read_text())
        raw = work / "raw"
        deposit = raw / (pdb + ".pdb")
        ligand = raw / (ccd + ".cif")
        lines = deposit.read_text().splitlines()
        row.update(
            pdb_url=f"https://files.rcsb.org/download/{pdb}.pdb",
            ccd_url=f"https://files.rcsb.org/ligands/download/{ccd}.cif",
            pdb_sha256=sha256(deposit),
            ccd_sha256=sha256(ligand),
            organism=next(
                l.split("ORGANISM_SCIENTIFIC: ")[1].split(";")[0]
                for l in lines
                if "ORGANISM_SCIENTIFIC:" in l
            ),
            doi=next(
                (
                    l.split("DOI")[1].strip().lower()
                    for l in lines
                    if l.startswith("JRNL        DOI")
                ),
                None,
            ),
        )
        if row["status"].startswith("Embedding/relaxation failed:"):
            row["status"] = (
                "Embedding or MMFF convergence failed for the single frozen conformer"
            )
        if pdb == "8K79":
            row["status"] += "; deposited chicken SRC"
        if row["status"] == "ready":
            row.update(crystal_contacts(gemmi.read_structure(str(deposit)), chain, ccd))
            row["bound_sdf_url"] = (
                f"https://models.rcsb.org/v1/{pdb.lower()}/ligand?auth_asym_id={chain}"
                f"&auth_seq_id={BOUND[pdb]}&encoding=sdf&copy_all_categories=false"
            )
            sdf = raw / (pdb + "-bound.sdf")
            if not sdf.exists():
                sdf.write_bytes(urlopen(row["bound_sdf_url"], timeout=60).read())
            row["bound_sdf_sha256"] = sha256(sdf)
            row["ligand_pdbqt_sha256"] = sha256(d / "ligand.pdbqt")
            row["receptor_pdbqt_sha256"] = sha256(d / "receptor.pdbqt")
        records.append(row)
    (work / "sources.json").write_text(json.dumps(records, indent=2) + "\n")
    (work / "software.json").write_text(
        json.dumps(
            {
                p: version(p)
                for p in [
                    "rdkit",
                    "meeko",
                    "gemmi",
                    "numpy",
                    "scipy",
                    "pdbfixer",
                    "openmm",
                ]
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("repo", type=Path)
    ap.add_argument("work", type=Path)
    a = ap.parse_args()
    main(a.repo.resolve(), a.work.resolve())
