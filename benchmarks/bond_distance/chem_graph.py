"""Check pair identifiers against deposited CCD bond tables or the input MOL2.

All graph distances use the chemical input. Vina supplies only the two serial numbers.
"""

import argparse
import csv
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

ORDERS = {
    "SING": Chem.BondType.SINGLE,
    "DOUB": Chem.BondType.DOUBLE,
    "TRIP": Chem.BondType.TRIPLE,
    "AROM": Chem.BondType.AROMATIC,
}
CASES = (
    ("5QCT", "E4Y", "tutorial/data/BACE_1_ligand.mol2"),
    ("1YET", "GDM", "hsp/raw/GDM.cif"),
    ("1FKF", "FK5", "extra/raw/FK5.cif"),
    ("1FKB", "RAP", "extra/raw/RAP.cif"),
)


def ccd(path):
    block = gemmi.cif.read(str(path)).sole_block()
    atoms = block.get_mmcif_category("_chem_comp_atom.")
    bonds = block.get_mmcif_category("_chem_comp_bond.")
    graph = Chem.RWMol()
    byname = {}
    ref_atom_stereo = {}
    ref_bond_stereo = {}
    for name, element, charge, stereo in zip(
        *(atoms[k] for k in ("atom_id", "type_symbol", "charge", "pdbx_stereo_config")),
        strict=True,
    ):
        atom = Chem.Atom(element.title())
        atom.SetFormalCharge(int(charge))
        atom.SetProp("source_name", name)
        byname[name] = graph.AddAtom(atom)
        if stereo in ("R", "S"):
            ref_atom_stereo[name] = stereo
    for a, b, order, stereo in zip(
        *(
            bonds[k]
            for k in ("atom_id_1", "atom_id_2", "value_order", "pdbx_stereo_config")
        ),
        strict=True,
    ):
        graph.AddBond(byname[a], byname[b], ORDERS[order])
        if stereo in ("E", "Z"):
            ref_bond_stereo[a, b] = stereo
    return Chem.RemoveHs(graph), ref_atom_stereo, ref_bond_stereo


def check(ligand, source):
    lines = ligand.read_text().splitlines()
    xyz, elements, closure = {}, {}, []
    for s in lines:
        if not s.startswith(("ATOM  ", "HETATM")):
            continue
        serial = int(s[6:11])
        ad = s.split()[-1]
        if ad.startswith(("G", "H")):
            continue
        element = {
            "A": "C",
            "C": "C",
            "N": "N",
            "NA": "N",
            "O": "O",
            "OA": "O",
            "S": "S",
            "SA": "S",
            "F": "F",
            "Cl": "Cl",
        }.get(ad)
        if ad.startswith("CG"):
            element = "C"
            closure.append(serial)
        assert element
        elements[serial] = element
        xyz[serial] = np.array([float(s[i : i + 8]) for i in (30, 38, 46)])
    if source.suffix == ".mol2":
        graph = Chem.MolFromMol2File(str(source), removeHs=True)
        for a in graph.GetAtoms():
            a.SetProp("source_name", a.GetProp("_TriposAtomName"))
        reference = Chem.Mol(graph)
        Chem.AssignStereochemistryFrom3D(reference)
        rdCIPLabeler.AssignCIPLabels(reference)
        atom_stereo = {
            a.GetProp("source_name"): a.GetProp("_CIPCode")
            for a in reference.GetAtoms()
            if a.HasProp("_CIPCode")
        }
        bond_stereo = {}
        refpos = graph.GetConformer().GetPositions()
        mapping = {}
        for j, p in xyz.items():
            (match,) = np.where(np.linalg.norm(refpos - p, axis=1) < 0.001)[0]
            mapping[j] = int(match)
    else:
        graph, atom_stereo, bond_stereo = ccd(source)
        smiles = Chem.MolFromSmiles(
            next(
                s.split()[2]
                for s in lines
                if s.startswith("REMARK SMILES ")
                and not s.startswith("REMARK SMILES IDX ")
            )
        )
        matches = graph.GetSubstructMatches(smiles, uniquify=False, useChirality=False)
        if len(matches) != 1:
            raise ValueError(
                f"Non-unique mapping of the complete chemical graph: {source}, {len(matches)}"
            )
        match = matches[0]
        numbers = [
            int(v)
            for s in lines
            if s.startswith("REMARK SMILES IDX ")
            for v in s.split()[3:]
        ]
        mapping = {
            s: match[i - 1] for i, s in zip(numbers[::2], numbers[1::2], strict=True)
        }
    assert set(mapping) == set(xyz) and set(mapping.values()) == set(
        range(graph.GetNumAtoms())
    )
    conf = Chem.Conformer(graph.GetNumAtoms())
    for serial, i in mapping.items():
        assert graph.GetAtomWithIdx(i).GetSymbol() == elements[serial]
        conf.SetAtomPosition(i, xyz[serial])
    graph.RemoveAllConformers()
    graph.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(graph)
    rdCIPLabeler.AssignCIPLabels(graph)
    atoms = {a.GetProp("source_name"): a for a in graph.GetAtoms()}
    mismatches = []
    for n, expected in atom_stereo.items():
        actual = atoms[n].GetProp("_CIPCode") if atoms[n].HasProp("_CIPCode") else "?"
        if actual != expected:
            mismatches.append(f"{n}:{expected}/{actual}")
    for (a, b), expected in bond_stereo.items():
        bond = graph.GetBondBetweenAtoms(atoms[a].GetIdx(), atoms[b].GetIdx())
        actual = bond.GetProp("_CIPCode") if bond.HasProp("_CIPCode") else "?"
        if actual != expected:
            mismatches.append(f"{a}-{b}:{expected}/{actual}")
    if mismatches:
        raise ValueError(
            f"Stereochemistry differs from source for {ligand}: {mismatches}"
        )
    end_names = {
        graph.GetAtomWithIdx(mapping[j]).GetProp("source_name") for j in closure
    }
    assert len(closure) == 2 and graph.GetBondBetweenAtoms(
        *(mapping[j] for j in closure)
    )
    return graph, mapping, end_names, len(atom_stereo), len(bond_stereo), mismatches


def main(work, pairfile, out):
    with pairfile.open() as f:
        records = list(csv.DictReader(f))
    with out.open("w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(
            (
                "pdb",
                "component",
                "source",
                "mapped_heavy",
                "stereocenters_checked",
                "double_bonds_checked",
                "stereo_mismatch_expected_observed",
                "serial_a",
                "serial_b",
                "source_path",
                "bonds",
                "cut_bond",
            )
        )
        if len(records) != 13:
            raise ValueError(
                "Expected the 13 parent pairs from the four frozen macrocycles"
            )
        for pdb, lig, relative in CASES:
            rows = [
                r
                for r in records
                if (pdb if pdb != "5QCT" else "BACE_1") in r["ligand"]
            ]
            path = work / rows[0]["ligand"]
            graph, mapping, ends, centers, alkenes, mismatches = check(
                path, work / relative
            )
            for row in rows:
                a, b = map(int, (row["serial_a"], row["serial_b"]))
                route = Chem.GetShortestPath(graph, mapping[a], mapping[b])
                names = [graph.GetAtomWithIdx(j).GetProp("source_name") for j in route]
                assert len(route) == 4 and any(
                    {i, j} == ends for i, j in zip(names, names[1:])
                )
                w.writerow(
                    (
                        pdb,
                        lig,
                        relative,
                        len(mapping),
                        centers,
                        alkenes,
                        ";".join(mismatches),
                        a,
                        b,
                        "-".join(names),
                        len(names) - 1,
                        "-".join(sorted(ends)),
                    )
                )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "parent_pairs", "out"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    main(a.work, a.parent_pairs, a.out)
