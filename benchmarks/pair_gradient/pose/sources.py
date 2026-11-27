"""Audit crystal contacts, input provenance and (optionally) shared preparation time."""

import argparse
import hashlib
import json
import shutil
import tempfile
import time
from pathlib import Path

import gemmi
import numpy as np
from paired import atoms, renin_correspondence, write
from prepare import CASES, base, prepare
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.spatial.distance import cdist


def crystal_contacts(structure, chain, component):
    model = structure[0]
    ns = gemmi.NeighborSearch(model, structure.cell, 8).populate(include_h=False)
    ligand = next(r for r in model[chain] if r.name == component)
    contacts = []
    for atom in ligand:
        if atom.element.name in ("H", "D"):
            continue
        for mark in ns.find_neighbors(atom, 0.01, 8):
            cra = mark.to_cra(model)
            if cra.residue.name == "HOH":
                continue
            pos = cra.atom.pos
            if mark.image_idx == 0 and atom.pos.dist(pos) <= 8:
                continue
            image = structure.cell.find_nearest_pbc_position(
                atom.pos, pos, mark.image_idx
            )
            contacts.append(
                (
                    atom.pos.dist(image),
                    atom.name,
                    cra.chain.name,
                    str(cra.residue.seqid),
                    cra.residue.name,
                    cra.atom.name,
                )
            )
    return {
        "mate_heavy_pairs_8A": len(contacts),
        "closest_mate_pair": min(contacts) if contacts else None,
    }


def reference_geometry(work):
    metadata = {
        p: json.loads((work / p / "preparation.json").read_text()) for p in CASES
    }
    renin_correspondence(work, metadata)
    rows = []
    for q, query in metadata.items():
        query_atoms = atoms((work / q / "ligand.pdbqt").read_text().splitlines())
        qrad = np.array([query_atoms[i][1] for i in query["serial"]])
        for r, rec in metadata.items():
            if rec["target"] != query["target"] or rec["receptor_status"] != "ready":
                continue
            rot, trans, ca_rmsd, count = base.alignment(query, rec)
            ref = np.array(query["mapped"][0]) @ rot + trans
            receptor = atoms((work / r / "receptor.pdbqt").read_text().splitlines())
            rxyz = np.array([v[0] for v in receptor.values()])
            rrad = np.array([v[1] for v in receptor.values()])
            d = cdist(ref, rxyz)
            rows.append(
                {
                    "query": q,
                    "receiver": r,
                    "kind": "cognate" if r == q else "cross",
                    "pocket_ca": count,
                    "pocket_rmsd_A": ca_rmsd,
                    "min_distance_A": d.min(),
                    "clashes": int(np.sum(d < 0.75 * (qrad[:, None] + rrad))),
                }
            )
    write(work / "reference_geometry.csv", rows)


def main(work, replay):
    reference_geometry(work)
    sources = []
    with tempfile.TemporaryDirectory() as temp:
        new = Path(temp)
        if replay:
            shutil.copytree(work / "raw", new / "raw")
        costs = []
        for pdb, (target, component, chain) in CASES.items():
            meta = json.loads((work / pdb / "preparation.json").read_text())
            structure = gemmi.read_structure(str(work / "raw" / (pdb + ".pdb")))
            mol = Chem.MolFromSmiles(meta["canonical_smiles"])
            files = [str(p.relative_to(work)) for p in (work / pdb).glob("*.pdbqt")]
            files += ["raw/" + pdb + ".pdb", "raw/" + component + ".cif"]
            digests = {
                f: hashlib.sha256((work / f).read_bytes()).hexdigest()
                for f in sorted(files)
            }
            sources.append(
                dict(
                    target=target,
                    pdb=pdb,
                    component=component,
                    chain=chain,
                    resolution_A=structure.resolution,
                    receptor_status=meta["receptor_status"],
                    nearest_asu_omitted=meta["closest_omitted_heavy_atom"],
                    nearest_rebuilt_A=meta.get("closest_repaired_ligand_distance"),
                    structure="https://www.rcsb.org/structure/" + pdb,
                    ccd="https://files.rcsb.org/ligands/download/" + component + ".cif",
                    smiles=meta["canonical_smiles"],
                    murcko=MurckoScaffold.MurckoScaffoldSmiles(mol=mol),
                    heavy_atoms=mol.GetNumHeavyAtoms(),
                    formal_charge=Chem.GetFormalCharge(mol),
                    **crystal_contacts(structure, chain, component),
                    sha256=digests,
                )
            )
            if replay:
                start = time.perf_counter()
                prepare(new, pdb)
                elapsed = time.perf_counter() - start
                matches = {
                    Path(f).name: (new / f).read_bytes() == (work / f).read_bytes()
                    for f in files
                    if f.endswith(".pdbqt")
                }
                if not all(matches.values()):
                    raise ValueError((pdb, matches))
                costs.append(
                    {"pdb": pdb, "seconds": elapsed, "byte_identical": matches}
                )
    (work / "sources.json").write_text(json.dumps(sources, indent=2) + "\n")
    if replay:
        (work / "preparation_replay.json").write_text(
            json.dumps(costs, indent=2) + "\n"
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--reprepare", action="store_true")
    a = p.parse_args()
    main(a.work, a.reprepare)
