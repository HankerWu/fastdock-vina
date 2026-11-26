"""Staurosporine whole-ligand RMSD, with only the receptor aligned in cross-docking."""

import argparse
import csv
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem

from freeze import SEEDS


def xyz(line):
    return [float(line[p : p + 8]) for p in (30, 38, 46)]


def reference(work):
    block = gemmi.cif.read_file(str(work / "raw/STU.cif")).sole_block()
    mol = Chem.RWMol()
    ids = {
        a: mol.AddAtom(Chem.Atom(e.title()))
        for a, e in block.find("_chem_comp_atom.", ["atom_id", "type_symbol"])
        if e != "H"
    }
    for a, b, order in block.find(
        "_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]
    ):
        if a in ids and b in ids:
            mol.AddBond(
                ids[a],
                ids[b],
                dict(SING=Chem.BondType.SINGLE, DOUB=Chem.BondType.DOUBLE)[order],
            )
    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    conf = Chem.Conformer(len(ids))
    for line in (work / "raw/1Q3D.pdb").read_text().splitlines():
        if line.startswith("HETATM") and line[17:22] == "STU A":
            conf.SetAtomPosition(ids[line[12:16].strip()], xyz(line))
    mol.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(mol)
    return mol


def evaluate(work, seeds=SEEDS):
    crystal = reference(work)
    template = (work / "prepared/44259.pdbqt").read_text().splitlines()
    mol = Chem.MolFromSmiles(
        next(
            l.split()[2]
            for l in template
            if l.startswith("REMARK SMILES ") and "IDX" not in l
        )
    )
    assert Chem.MolToInchiKey(mol) == Chem.MolToInchiKey(crystal)
    maps = crystal.GetSubstructMatches(mol, uniquify=False, useChirality=True)
    assert maps and mol.GetNumAtoms() == 35
    ids = {}
    for line in template:
        if line.startswith("REMARK SMILES IDX"):
            values = list(map(int, line.split()[3:]))
            ids.update({i - 1: s for i, s in zip(values[::2], values[1::2])})
    xyz_ref = crystal.GetConformer().GetPositions()

    def cas(code):
        return {
            line[22:27]: np.array(xyz(line))
            for line in (work / f"raw/{code}.pdb").read_text().splitlines()
            if line.startswith("ATOM  ")
            and line[12:16].strip() == "CA"
            and line[21] == "A"
            and line[16] in (" ", "A")
        }

    rows = []
    for receptor in ("1Q3D", "1Q5K"):
        ca = cas("1Q3D")
        cb = cas(receptor)
        common = sorted(
            i
            for i in ca.keys() & cb
            if np.min(np.linalg.norm(ca[i] - xyz_ref, axis=1)) <= 10
        )
        a, b = [np.array([c[i] for i in common]) for c in (ca, cb)]
        u, _, vh = np.linalg.svd((a - a.mean(0)).T @ (b - b.mean(0)))
        rot = u @ np.diag([1, 1, np.sign(np.linalg.det(u @ vh))]) @ vh
        translate = b.mean(0) - a.mean(0) @ rot

        def rms(x, y):
            return float(np.sqrt(np.sum((x - y) ** 2) / len(x)))

        pocket_rmsd = rms(a @ rot + translate, b)
        target = xyz_ref @ rot + translate
        for arm, folder in [
            ("native", "pools"),
            ("mixed", "mixed-pools"),
            ("soft", "soft-pools"),
        ]:
            for seed in seeds:
                path = work / f"{folder}/{receptor}_44259_{seed}.pdbqt"
                for model in path.read_text().split("ENDMDL"):
                    if not model.strip():
                        continue
                    lines = model.strip().splitlines()
                    coordinates = {
                        int(l[6:11]): xyz(l)
                        for l in lines
                        if l.startswith(("ATOM  ", "HETATM"))
                    }
                    ligand = np.array(
                        [coordinates[ids[i]] for i in range(mol.GetNumAtoms())]
                    )
                    rows.append(
                        dict(
                            arm=arm,
                            receptor=receptor,
                            seed=seed,
                            model=int(lines[0].split()[1]),
                            heavy=35,
                            mappings=len(maps),
                            aligned_ca=len(a),
                            pocket_rmsd=pocket_rmsd,
                            rmsd=min(rms(ligand, target[list(m)]) for m in maps),
                        )
                    )
    with (work / "results/stu-rmsd.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    args = p.parse_args()
    evaluate(args.work, args.seeds)
