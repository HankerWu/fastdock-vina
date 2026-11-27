"""Prepare the two independent human VHL structures and CCD-graph ligands."""

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import gemmi
import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import cdist

BENCH = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(BENCH / "paired_crossdocking"))
import run as paired

CASES = {
    "4W9H": (
        "3JF",
        "Cc1c(scn1)c2ccc(cc2)CNC(=O)[C@@H]3C[C@H](CN3C(=O)[C@H](C(C)(C)C)NC(=O)C)O",
        "GFVIEZBZIUKYOG-SVFBPWRDSA-N",
    ),
    "5LLI": (
        "6Z3",
        "Cc1c(scn1)c2ccc(cc2)CNC(=O)[C@@H]3C[C@H](CN3C(=O)[C@H](C(C)(C)C)NC(=O)C4(CC4)C#N)O",
        "NDVQUNZCNAMROD-RZUBCFFCSA-N",
    ),
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def xyz(line):
    return [float(line[i : i + 8]) for i in (30, 38, 46)]


def ligand(folder, lines, atoms, smiles, key):
    mol = Chem.MolFromSmiles(smiles)
    assert Chem.MolToInchiKey(mol) == key
    h = Chem.AddHs(mol)
    settings = AllChem.ETKDGv3()
    settings.randomSeed = 51127
    assert AllChem.EmbedMolecule(h, settings) == 0
    assert AllChem.MMFFOptimizeMolecule(h, maxIters=2000) == 0
    (setup,) = MoleculePreparation().prepare(h)
    txt, ok, error = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(error)
    (folder / "ligand.pdbqt").write_text(txt)
    serial = {}
    for line in txt.splitlines():
        if line.startswith("REMARK INDEX MAP"):
            ids = list(map(int, line.split()[3:]))
            serial.update(zip(ids[::2], ids[1::2]))
    serial = [serial[i + 1] for i in range(mol.GetNumAtoms())]

    ids = {int(l[6:11]) for l in atoms}
    bonds = [l for l in lines if l.startswith("CONECT") and int(l[6:11]) in ids]
    block = "\n".join([*atoms, *bonds, "END"]) + "\n"
    raw = Chem.MolFromPDBBlock(
        block, sanitize=False, removeHs=False, proximityBonding=False
    )
    ref = AllChem.AssignBondOrdersFromTemplate(mol, raw)
    Chem.AssignStereochemistryFrom3D(ref, replaceExistingTags=True)
    assert Chem.MolToInchiKey(ref) == key and ref.GetNumAtoms() == mol.GetNumAtoms()
    matches = ref.GetSubstructMatches(
        mol, uniquify=False, useChirality=True, maxMatches=10001
    )
    assert 0 < len(matches) < 10001
    return serial, [
        ref.GetConformer().GetPositions()[list(m)].tolist() for m in matches
    ]


def prepare(work, pdb):
    code, smiles, key = CASES[pdb]
    path = work / "raw" / (pdb + ".pdb")
    lines = path.read_text().splitlines()
    atoms = [
        l
        for l in lines
        if l.startswith(("ATOM  ", "HETATM"))
        and l[76:78].strip() not in ("H", "D")
        and l[16] in " A"
    ]
    ref_atoms = [l for l in atoms if l[17:20] == code and l[21:26] == "C 301"]
    folder = work / pdb
    folder.mkdir(exist_ok=True)
    serial, refs = ligand(folder, lines, ref_atoms, smiles, key)
    crystal = np.array(refs[0])

    def nearest(records):
        return (
            float(cdist(crystal, [xyz(l) for l in records]).min()) if records else None
        )

    groups = {
        "VHL_CAS": [l for l in atoms if l[17:20] == "CAS" and l[21] == "C"],
        "Elongins": [l for l in atoms if l.startswith("ATOM") and l[21] in "AB"],
        "water": [l for l in atoms if l[17:20] == "HOH"],
        "other_asu_protein": [
            l for l in atoms if l.startswith("ATOM") and l[21] not in "ABC"
        ],
        "other_nonwater": [
            l
            for l in atoms
            if l.startswith("HETATM") and l[17:20] not in [code, "HOH", "CAS"]
        ],
        "other_ligand_copy": [l for l in atoms if l[17:20] == code and l[21] != "C"],
    }
    distances = {k: nearest(v) for k, v in groups.items()}
    assert distances["other_nonwater"] is None or distances["other_nonwater"] > 4
    protein = []
    for l in atoms:
        if l[21] not in "ABC":
            continue
        if l.startswith("ATOM"):
            protein.append(l)
        elif l[17:20] == "CAS" and l[12:16].strip() in (
            "N",
            "CA",
            "C",
            "O",
            "CB",
            "SG",
        ):
            protein.append("ATOM  " + l[6:17] + "CYS" + l[20:])
    # Both side chains on 4W9H C142/C143 end at CB in the deposit. PDBFixer
    # rebuilds mutually overlapping atoms which Meeko connects covalently.
    # Keep their deposited heavy atoms, using the matching alanine template.
    truncations = {}
    for index in (142, 143):
        records = [l for l in protein if l[21] == "C" and int(l[22:26]) == index]
        if {l[12:16].strip() for l in records} == {"N", "CA", "C", "O", "CB"}:
            assert nearest(records) > 12
            truncations[str(index)] = [records[0][17:20], nearest(records)]
            protein = [l[:17] + "ALA" + l[20:] if l in records else l for l in protein]
    all_atoms, repaired = [], {}
    for ch in "ABC":
        single = folder / ch
        single.mkdir(exist_ok=True)
        chosen = [l[:21] + "A" + l[22:] for l in protein if l[21] == ch]
        repaired[ch] = paired.receptor(chosen, single, pdb, crystal)
        for l in (single / "receptor.pdbqt").read_text().splitlines():
            if l.startswith(("ATOM", "HETATM")):
                all_atoms.append(
                    l[:6] + f"{len(all_atoms) + 1:5d}" + l[11:21] + ch + l[22:]
                )
    (folder / "receptor.pdbqt").write_text("\n".join(all_atoms) + "\n")
    ca = {
        l[21:27] + l[17:20]: xyz(l)
        for l in protein
        if l[12:16].strip() == "CA" and nearest([l]) <= 12
    }
    # Assess lattice contacts including periodic translations outside the ASU.
    sys.path.insert(0, str(BENCH / "pair_gradient/pose"))
    from sources import crystal_contacts

    struct = gemmi.read_structure(str(path))
    m = dict(
        pdb=pdb,
        component=code,
        smiles=smiles,
        inchikey=key,
        crystal=refs,
        serial=serial,
        center=crystal.mean(axis=0).tolist(),
        ca=ca,
        closest_A=distances,
        truncated_at_deposited_CB=truncations,
        repair=repaired,
        min_ligand_occupancy=min(float(l[54:60]) for l in ref_atoms),
        resolution_A=struct.resolution,
        raw_sha256=digest(path),
        receptor_sha256=digest(folder / "receptor.pdbqt"),
        ligand_sha256=digest(folder / "ligand.pdbqt"),
        **crystal_contacts(struct, "C", code),
    )
    (folder / "metadata.json").write_text(json.dumps(m, indent=2) + "\n")
    print(pdb, len(serial), len(refs), distances, flush=True)


if __name__ == "__main__":
    app = argparse.ArgumentParser(description=__doc__)
    app.add_argument("work", type=Path)
    args = app.parse_args()
    os.environ["OPENMM_CPU_THREADS"] = "1"
    for key in CASES:
        prepare(args.work, key)
