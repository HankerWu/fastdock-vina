"""Cross-dock the one-sided 1FAX ligand into two independently prepared FXa receptors."""

import argparse
import itertools
from pathlib import Path
import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from experiment import ARMS, SEEDS, prepare, run, write
from contrasts import EXTRA
from evaluate import pose_models
from make_ion_pair_anchors import carboxylate, ligand_sites, xyz
from score_physical_pool import serial_index_pairs


def structure_lines(path):
    return [
        s
        for s in path.read_text().splitlines()
        if s.startswith(("ATOM  ", "HETATM")) and s[21] == "A" and s[16] in " A"
    ]


def chemical_graph(pdb, cif):
    block = gemmi.cif.read_file(str(cif)).sole_block()
    row = block.find("_chem_comp_atom.", ["atom_id", "type_symbol"])
    atoms = {name: element for name, element in row if element != "H"}
    points = {
        l[12:16].strip(): xyz(l)
        for l in structure_lines(pdb)
        if l[17:20] == "DX9" and l[76:78].strip() != "H"
    }
    assert atoms.keys() == points.keys()
    mol = Chem.RWMol()
    ids = {
        name: mol.AddAtom(Chem.Atom(element.title())) for name, element in atoms.items()
    }
    for a, b, order in block.find(
        "_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]
    ):
        if a in ids and b in ids:
            mol.AddBond(
                ids[a],
                ids[b],
                {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE}[order],
            )
    lig = mol.GetMol()
    Chem.SanitizeMol(lig)
    conf = Chem.Conformer(len(ids))
    for n, i in ids.items():
        conf.SetAtomPosition(i, points[n])
    lig.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(lig)
    neutral = Chem.MolToSmiles(lig)
    (canonical,) = [
        d.strip('"')
        for kind, software, d in block.find(
            "_pdbx_chem_comp_descriptor.", ["type", "program", "descriptor"]
        )
        if kind == "SMILES_CANONICAL" and software == "CACTVS"
    ]
    if neutral != Chem.MolToSmiles(Chem.MolFromSmiles(canonical)):
        raise ValueError(
            "Ligand stereochemistry disagrees between the CCD and the crystal instance"
        )
    for atom in lig.GetAtoms():
        if (
            atom.GetAtomicNum() == 7
            and atom.GetDegree() == 1
            and atom.GetBonds()[0].GetBondType() == Chem.BondType.DOUBLE
        ):
            atom.SetFormalCharge(1)
        if (
            atom.GetAtomicNum() == 8
            and atom.GetDegree() == 1
            and atom.GetBonds()[0].GetBondType() == Chem.BondType.SINGLE
        ):
            atom.SetFormalCharge(-1)
    Chem.SanitizeMol(lig)
    assert Chem.GetFormalCharge(lig) == 1
    return lig, np.array([points[n] for n in ids]), neutral


def alignment(query, target, lasp):
    def ca(path):
        return {
            (s[22:26], s[26]): np.array(xyz(s))
            for s in structure_lines(path)
            if s.startswith("ATOM  ") and s[12:16].strip() == "CA"
        }

    x, y = ca(query), ca(target)
    shared = [i for i in x.keys() & y.keys() if np.linalg.norm(x[i] - lasp) < 15]
    shared.sort()
    xx = np.array([x[i] for i in shared])
    yy = np.array([y[i] for i in shared])
    u, _, v = np.linalg.svd((xx - xx.mean(0)).T @ (yy - yy.mean(0)))
    scaling = np.diag([1, 1, np.sign(np.linalg.det(u @ v))])
    rot = u @ scaling @ v
    shift = yy.mean(0) - xx.mean(0) @ rot
    return (
        rot,
        shift,
        len(shared),
        float(np.sqrt(np.square(xx @ rot + shift - yy).sum() / len(xx))),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "fjs", "p16", "binary", "result"):
        p.add_argument(name, type=Path)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--collect", action="store_true")
    a = p.parse_args()
    raw = a.work / "raw"
    out = a.work / "cross"
    out.mkdir(exist_ok=True)
    crystal = raw / "1FAX.pdb"
    graph, ref, neutral = chemical_graph(crystal, raw / "DX9.cif")
    # Release all crystal coordinates before preparing the input by graph alone.
    graph.RemoveAllConformers()
    smiles = Chem.MolToSmiles(graph)
    print(
        "1FAX coordinate-verified neutral and docked graph", neutral, smiles, flush=True
    )
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    par = AllChem.ETKDGv3()
    par.randomSeed = 2026
    if AllChem.EmbedMolecule(m, par) or AllChem.MMFFOptimizeMolecule(m, maxIters=1000):
        raise ValueError("Exposure of an unconverged structure")
    (setup,) = MoleculePreparation().prepare(m)
    text, ok, why = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(why)
    ligand = out / "DX9.pdbqt"
    ligand.write_text(text)
    inputs = {}
    sources = {"1FJS": a.fjs, "2P16": a.p16}
    for rec, folder in sources.items():
        prepare(out, rec, folder / f"{rec}.pdb", "A", [ligand])
        inputs[rec, "DX9"] = (ligand, folder / f"{rec}.pdbqt", folder / "box.txt")
    cases = list(itertools.product(sources, ("DX9",), (*SEEDS, *EXTRA), ARMS))
    if not a.collect:
        run(a.binary.resolve(), out, inputs, cases, a.jobs)
    maps = dict(
        (i - 1, s)
        for i, s in serial_index_pairs(text.splitlines(), "REMARK SMILES IDX ")
    )
    smi = next(
        l.split()[2]
        for l in text.splitlines()
        if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX")
    )

    # All true symmetry permutations, including resonance-carboxylate and amidinium swaps.
    def resonance_copy(mol):
        mol = Chem.RWMol(mol)
        for q in ("[C](=[N+;H2])-[N;H2]", "[C](=O)-[O-]"):
            for c, a, b in mol.GetSubstructMatches(Chem.MolFromSmarts(q)):
                for i in (a, b):
                    m = mol.GetAtomWithIdx(i)
                    m.SetFormalCharge(0)
                    m.SetNumExplicitHs(0)
                    m.SetNoImplicit(True)
                    mol.GetBondBetweenAtoms(c, i).SetBondType(Chem.BondType.SINGLE)
        return mol.GetMol()

    order = resonance_copy(graph).GetSubstructMatches(
        resonance_copy(Chem.MolFromSmiles(smi)), uniquify=False, useChirality=False
    )
    assert len(order) == 8, len(order)
    lasp = np.array(
        next(
            xyz(l)
            for l in structure_lines(crystal)
            if l[12:16].strip() == "CG" and l[22:26].strip() == "189"
        )
    )
    rows = []
    ((_, n1, n2),) = ligand_sites(ligand)[0]
    for rec, folder in sources.items():
        rot, shift, count, err = alignment(crystal, folder / f"{rec}.pdb", lasp)
        target = ref @ rot + shift
        oxy = carboxylate(folder / f"{rec}.pdbqt", "A", 189)[1:]
        for seed, arm in itertools.product((*SEEDS, *EXTRA), ARMS):
            rms = []
            contacts = []
            for _, x in pose_models(out / f"{rec}_DX9_{seed}_{arm}.pdbqt"):
                pred = np.array([x[maps[i]] for i in range(len(maps))])
                contact = min(
                    (
                        sorted(np.linalg.norm(np.array([x[n1], x[n2]]) - os, axis=1))
                        for os in (oxy, oxy[::-1])
                    ),
                    key=lambda ds: max(ds),
                )
                contacts.append(contact)
                rms.append(
                    min(
                        float(
                            np.sqrt(np.sum((pred - target[list(o)]) ** 2) / len(pred))
                        )
                        for o in order
                    )
                )
            rows.append(
                dict(
                    receptor=rec,
                    seed=seed,
                    arm=arm,
                    CA_pairs=count,
                    CA_rms=err,
                    rms_top=rms[0],
                    rms_best=min(rms),
                    n_modes=len(rms),
                    top_short=contacts[0][0],
                    top_long=contacts[0][1],
                    oracle_long=contacts[int(np.argmin(rms))][1],
                )
            )
    write(a.result, rows)


if __name__ == "__main__":
    main()
