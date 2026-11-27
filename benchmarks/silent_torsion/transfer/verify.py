"""Rebuild reference graphs and pocket transforms from depositions and audit every RMSD."""

import argparse
import csv
import itertools
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SVDSuperimposer import SVDSuperimposer
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdMolAlign

COMPONENTS = {"2LBD": "REA", "1FBY": "9CR", "1Q41": "IXM", "1Q5K": "TMU"}


def deposition(work, pdb, code):
    chain = PDBParser(QUIET=True).get_structure(pdb, work / "raw" / (pdb + ".pdb"))[0][
        "A"
    ]
    (ligand,) = (r for r in chain if r.get_resname() == code)
    cif = MMCIF2Dict(work / "raw" / (code + ".cif"))
    g = Chem.RWMol()
    atom_names = {}
    for name, element, charge in zip(
        *(cif["_chem_comp_atom." + k] for k in ("atom_id", "type_symbol", "charge"))
    ):
        if element.upper() in {"H", "D"}:
            continue
        atom = Chem.Atom(element.title())
        atom.SetFormalCharge(int(charge))
        atom_names[name] = g.AddAtom(atom)
    for n, m, order in zip(
        *(
            cif["_chem_comp_bond." + k]
            for k in ("atom_id_1", "atom_id_2", "value_order")
        )
    ):
        if n in atom_names and m in atom_names:
            g.AddBond(
                atom_names[n],
                atom_names[m],
                {
                    "SING": Chem.BondType.SINGLE,
                    "DOUB": Chem.BondType.DOUBLE,
                    "TRIP": Chem.BondType.TRIPLE,
                    "AROM": Chem.BondType.AROMATIC,
                }[order],
            )
    geom = g.GetMol()
    Chem.SanitizeMol(geom)
    ligand_atoms = {a.get_name(): a.get_coord().astype(float) for a in ligand}
    if set(atom_names) != set(ligand_atoms) or any(a.is_disordered() for a in ligand):
        raise ValueError(
            f"{pdb}: incomplete or alternate-position reference requires review"
        )
    conf = Chem.Conformer(geom.GetNumAtoms())
    for name, i in atom_names.items():
        conf.SetAtomPosition(i, ligand_atoms[name])
    geom.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(geom)
    center = np.mean(list(ligand_atoms.values()), axis=0)
    backbone = {
        (r.id, r.get_resname()): r["CA"].get_coord().astype(float)
        for r in chain
        if r.id[0] == " "
        and "CA" in r
        and np.linalg.norm(r["CA"].get_coord() - center) < 14
    }
    return geom, backbone


def transformed_reference(source, receiver):
    lig, ca = source
    _, target_ca = receiver
    keys = sorted(ca.keys() & target_ca.keys())
    align = SVDSuperimposer()
    align.set(np.array([target_ca[k] for k in keys]), np.array([ca[k] for k in keys]))
    align.run()
    rotation, translation = align.get_rotran()
    ref = Chem.Mol(lig)
    conf = ref.GetConformer()
    for i, xyz in enumerate(conf.GetPositions() @ rotation + translation):
        conf.SetAtomPosition(i, xyz)
    return ref


def main(work):
    data = {p: deposition(work, p, c) for p, c in COMPONENTS.items()}
    frame = {
        (q, r): transformed_reference(data[q], data[r])
        for q, r in itertools.product(data, repeat=2)
        if q == r or {q, r} == {"1Q41", "1Q5K"}
    }
    check = {}
    with (work / "modes.csv").open() as f:
        for row in csv.DictReader(f):
            check.setdefault(row["key"], {})[int(row["mode"])] = float(row["rmsd_A"])
    maximum, total = 0.0, 0
    for key, rms in check.items():
        q, r = key.split("_")[:2]
        ref = frame[q, r]
        pdbqt = PDBQTMolecule.from_file(
            str(work / "runs" / (key + ".pdbqt")), skip_typing=True
        )
        (mol,) = RDKitMolCreate.from_pdbqt_mol(pdbqt)
        probe = Chem.RemoveHs(mol)
        matches = ref.GetSubstructMatches(
            probe, uniquify=False, useChirality=True, maxMatches=10001
        )
        if not matches or len(matches) > 10000 or probe.GetNumConformers() != len(rms):
            raise ValueError(key)
        mappings = [list(enumerate(m)) for m in matches]
        for i, conf in enumerate(probe.GetConformers(), 1):
            value = rdMolAlign.CalcRMS(
                probe,
                ref,
                prbId=conf.GetId(),
                map=mappings,
                symmetrizeConjugatedTerminalGroups=False,
            )
            maximum = max(maximum, abs(value - rms[i]))
            total += 1
    if not total or maximum > 1e-4:
        raise ValueError((total, maximum))
    print(
        f"{len(check)} pools, {total} modes, maximum absolute RMSD difference {maximum:.9g} A"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    main(parser.parse_args().work)
