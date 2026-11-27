"""Reconstruct pose RMSDs from PDB/CCD and serialized output without preparation metadata."""

import argparse
import csv
import json
from pathlib import Path

import gemmi
import numpy as np
from panel import PANELS
from rdkit import Chem
from scipy.spatial.transform import Rotation


def reference(work, pdb, component, chain_id="A"):
    # Use the mmCIF dictionary bond table and PDB-deposited coordinates.
    chain = gemmi.read_structure(str(work / "raw" / (pdb + ".pdb")))[0][chain_id]
    ligand = next(r for r in chain if r.name == component)
    positions = {
        a.name: list(a.pos)
        for a in ligand
        if a.element.name not in ("H", "D") and a.altloc in ("\x00", " ", "A")
    }
    center = np.mean(list(positions.values()), axis=0)
    atoms = {}
    for residue in chain:
        ca = residue.find_atom("CA", "*")
        if ca is not None and residue.het_flag == "A":
            atoms[(residue.seqid.num, residue.seqid.icode, residue.name)] = np.array(
                list(ca.pos)
            )
    ccd = gemmi.cif.read_file(str(work / "raw" / (component + ".cif"))).sole_block()
    mol = Chem.RWMol()
    idx = {}
    for key, symbol, charge in ccd.find(
        "_chem_comp_atom.", ["atom_id", "type_symbol", "charge"]
    ):
        name = gemmi.cif.as_string(key)
        if symbol in ("H", "D"):
            continue
        atom = Chem.Atom(symbol.title())
        atom.SetFormalCharge(int(charge))
        idx[name] = mol.AddAtom(atom)
    for a, b, bond in ccd.find(
        "_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]
    ):
        a = gemmi.cif.as_string(a)
        b = gemmi.cif.as_string(b)
        if a in idx and b in idx:
            mol.AddBond(
                idx[a],
                idx[b],
                {
                    "SING": Chem.BondType.SINGLE,
                    "DOUB": Chem.BondType.DOUBLE,
                    "TRIP": Chem.BondType.TRIPLE,
                }[bond],
            )
    geom = mol.GetMol()
    Chem.SanitizeMol(geom)
    conf = Chem.Conformer(geom.GetNumAtoms())
    for name, i in idx.items():
        conf.SetAtomPosition(i, positions[name])
    geom.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(geom, replaceExistingTags=True)
    return atoms, center, geom


def verify(work):
    dic = {
        p: reference(work, p, c) for group in PANELS.values() for p, c in group.items()
    }
    archived = {}
    for arm in ("native", "cubic"):
        with (work / f"modes_{arm}_24000.csv").open() as f:
            for row in csv.DictReader(f):
                archived[
                    (row["query"], row["receiver"], row["seed"], arm, int(row["model"]))
                ] = float(row["rmsd"])
    deviations = []
    for file in sorted((work / "runs").glob("*.pdbqt")):
        q, r, seed, arm, _ = file.stem.split("_")
        qca, qcenter, mol = dic[q]
        rca, rcenter, _ = dic[r]
        keys = sorted(
            k
            for k in qca.keys() & rca.keys()
            if np.linalg.norm(qca[k] - qcenter) < 14
            and np.linalg.norm(rca[k] - rcenter) < 14
        )
        x = np.stack([qca[k] for k in keys])
        y = np.stack([rca[k] for k in keys])
        xm = x.mean(0)
        ym = y.mean(0)
        rotation, _ = Rotation.align_vectors(y - ym, x - xm)
        ref = rotation.apply(mol.GetConformer().GetPositions() - xm) + ym
        for num, block in enumerate(file.read_text().split("ENDMDL")[:-1], 1):
            lines = block.splitlines()
            smi = next(
                l.split(maxsplit=2)[2]
                for l in lines
                if l.startswith("REMARK SMILES ")
                and not l.startswith("REMARK SMILES IDX")
            )
            query = Chem.MolFromSmiles(smi)
            serial = {}
            for l in lines:
                if l.startswith("REMARK SMILES IDX"):
                    values = list(map(int, l.split()[3:]))
                    serial.update(zip(values[::2], values[1::2], strict=True))
            coordinates = {
                int(l[6:11]): [float(l[i : i + 8]) for i in (30, 38, 46)]
                for l in lines
                if l.startswith(("ATOM  ", "HETATM"))
            }
            coords = np.array(
                [coordinates[serial[i + 1]] for i in range(query.GetNumAtoms())]
            )
            maps = mol.GetSubstructMatches(
                query, uniquify=False, useChirality=True, maxMatches=10001
            )
            assert 0 < len(maps) < 10001
            rmsd = min(
                np.linalg.norm(ref[list(m)] - coords) / np.sqrt(len(m)) for m in maps
            )
            key = (q, r, seed, arm, num)
            deviations.append(abs(archived.pop(key) - rmsd))
    if archived or max(deviations) > 1e-5:
        raise ValueError((len(archived), max(deviations)))
    return {
        "verified_modes": len(deviations),
        "maximum_rmsd_difference_A": max(deviations),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    a = p.parse_args()
    result = verify(a.work)
    print(json.dumps(result))
    (a.work / "independent.json").write_text(json.dumps(result, indent=2) + "\n")
