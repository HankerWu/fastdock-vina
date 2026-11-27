"""Frozen panel and input preparation; common historical mapping and PDBQT recipe."""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import gemmi
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH / "pair_gradient/pose"))
import paired
import prepare
from sources import crystal_contacts

__all__ = ["ARMS", "CAPS", "CASES", "SEEDS", "paired"]

CASES = {
    "1K7L": ("PPARA", "544", "A"),
    "2P54": ("PPARA", "735", "A"),
    "4UYE": ("BRPF1", "9F9", "A"),
    "5FG4": ("BRPF1", "5XE", "A"),
    "1A28": ("PGR", "STR", "A"),
    "1SQN": ("PGR", "NDR", "A"),
}
CAPS = (3000, 12000)
SEEDS = (297131, 297149)
ARMS = ("upstream", "parent", "halving", "quadratic")


def prepare_input(work, pdb):
    prepare.CASES = CASES
    start = time.perf_counter()
    prepare.prepare(work, pdb)
    target, component, chain = CASES[pdb]
    meta = json.loads((work / pdb / "preparation.json").read_text())
    structure = gemmi.read_structure(str(work / "raw" / (pdb + ".pdb")))
    mol = Chem.MolFromSmiles(meta["canonical_smiles"])
    source = dict(
        target=target,
        pdb=pdb,
        ccd=component,
        chain=chain,
        structure_source=f"https://www.rcsb.org/structure/{pdb}",
        ccd_source=f"https://files.rcsb.org/ligands/download/{component}.cif",
        resolution_A=structure.resolution,
        formal_charge=Chem.GetFormalCharge(mol),
        murcko_smiles=MurckoScaffold.MurckoScaffoldSmiles(mol=mol),
        heavy_atoms=mol.GetNumHeavyAtoms(),
        seconds=time.perf_counter() - start,
        **crystal_contacts(structure, chain, component),
    )
    (work / pdb / "source.json").write_text(json.dumps(source, indent=2) + "\n")


def write_sources(work):
    rows = []
    for pdb in CASES:
        folder = work / pdb
        base = json.loads((folder / "preparation.json").read_text())
        source = json.loads((folder / "source.json").read_text())
        torsions = next(
            int(line.split()[1])
            for line in (folder / "ligand.pdbqt").read_text().splitlines()
            if line.startswith("TORSDOF")
        )
        rows.append(
            {
                **source,
                "status": base["receptor_status"],
                "nearest_repaired_A": base.get("closest_repaired_ligand_distance"),
                "nearest_omitted_nonwater": base["closest_omitted_heavy_atom"],
                "torsdof": torsions,
                "smiles": base["canonical_smiles"],
                "pdb_sha256": base["pdb_sha256"],
                "ccd_sha256": base["ccd_sha256"],
            }
        )
    (work / "sources.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument("ids", nargs="*")
    a = ap.parse_args()
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    for code in a.ids or CASES:
        prepare_input(a.work, code)
