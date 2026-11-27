"""Preserve deposited nearby buffer compounds in a separately labelled sensitivity model."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
from collections import defaultdict

import gemmi
import numpy as np
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

from prepare import panel, prepare_pose


def retained(root, pdb, source):
    lines = (root / (pdb + ".pdb")).read_text().splitlines()
    species = defaultdict(list)
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and line[76:78].strip() not in (
            "D",
            "H",
        ):
            species[line[17:27]].append(line)
    lig = species[source["ligand"].rjust(3) + " " + source["ligand_instance"]]
    observed = np.array([panel.base.coords(l) for l in lig])
    keep = {}
    for key, atoms in species.items():
        if (
            atoms == lig
            or key[:3] in ("HOH", "DOD", "WAT")
            or (atoms[0].startswith("ATOM") and atoms[0][21] == source["chain"])
        ):
            continue
        mindist = min(
            np.linalg.norm(observed - np.array(panel.base.coords(a)), axis=1).min()
            for a in atoms
        )
        if mindist < 8:
            if atoms[0].startswith("ATOM"):
                raise ValueError(f"another protein chain near ligand ({key})")
            for l in lines:
                if l.startswith("LINK") and key in (l[17:27], l[47:57]):
                    raise ValueError(f"deposited covalent LINK to {key}")
            keep[key] = (float(mindist), atoms)
    return lines, keep


def parametrize(ccd, atoms, idx):
    table = gemmi.cif.read_file(str(ccd)).sole_block()
    rec = table.get_mmcif_category("_chem_comp_atom.")
    graph = Chem.RWMol()
    atom_indices = {}
    for n, e, q in zip(rec["atom_id"], rec["type_symbol"], rec["charge"]):
        atom = Chem.Atom(e.title())
        atom.SetFormalCharge(int(q))
        atom.SetProp("ccd", n)
        atom_indices[n] = graph.AddAtom(atom)
    bonds = table.get_mmcif_category("_chem_comp_bond.")
    for a, b, t in zip(bonds["atom_id_1"], bonds["atom_id_2"], bonds["value_order"]):
        graph.AddBond(
            atom_indices[a],
            atom_indices[b],
            dict(SING=Chem.BondType.SINGLE, DOUB=Chem.BondType.DOUBLE)[t],
        )
    mol = graph.GetMol()
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    alt = sorted({l[16] for l in atoms} - {" "})
    lookup = {
        l[12:16].strip(): l for l in atoms if l[16] in (" ", alt[0] if alt else " ")
    }
    names = [a.GetProp("ccd") for a in mol.GetAtoms()]
    if set(names) != set(lookup):
        raise ValueError(f"incomplete additive {ccd.name}")
    conf = Chem.Conformer(len(names))
    for i, n in enumerate(names):
        conf.SetAtomPosition(i, panel.base.coords(lookup[n]))
    mol.AddConformer(conf)
    mol = Chem.AddHs(mol, addCoords=True)
    (setup,) = MoleculePreparation().prepare(mol)
    text, ok, error = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(error)
    mp = {}
    for l in text.splitlines():
        if l.startswith("REMARK INDEX MAP"):
            m = [int(x) for x in l.split()[3:]]
            mp.update(zip(m[1::2], (i - 1 for i in m[::2])))
    output = []
    for l in text.splitlines():
        if not l.startswith(("ATOM  ", "HETATM")):
            continue
        i = mp[int(l[6:11])]
        canonical = names[i] if i < len(names) else l[12:16].strip()
        if (
            i < len(names)
            and np.linalg.norm(
                np.asarray(panel.base.coords(l)) - panel.base.coords(lookup[canonical])
            )
            > 0.001
        ):
            raise ValueError("additive coordinate changed")
        source = atoms[0]
        output.append(
            "HETATM"
            + f"{idx + len(output):5d} "
            + canonical.rjust(4)
            + " "
            + source[17:27]
            + l[27:]
        )
    return output, Chem.GetFormalCharge(mol)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("prepared", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    src = a.prepared.resolve()
    w = a.out.resolve()
    w.mkdir(parents=True, exist_ok=True)
    original = json.loads((src / "manifest.json").read_text())
    results = []
    for row in original:
        if row["status"] != "excluded":
            continue
        pdb, code, target = (row[k] for k in ("pdb", "ligand", "target"))
        d = w / pdb
        d.mkdir(exist_ok=True)
        source = json.loads((src / pdb / "source.json").read_text())
        row = dict(row)
        try:
            lines, keep = retained(src / pdb, pdb, source)
            addons, detail = [], []
            for key, (distance, atoms) in keep.items():
                chem = key[:3].strip()
                ccd = d / (chem + ".cif")
                sha = panel.base.fetch(
                    ccd, f"https://files.rcsb.org/ligands/download/{chem}.cif"
                )
                atom_lines, charge = parametrize(ccd, atoms, 8000 + len(addons))
                addons.extend(atom_lines)
                detail.append(
                    dict(residue=key, distance=distance, charge=charge, sha256=sha)
                )
            # The general preparer constructs the protein-only part and checks every remaining omitted residue.
            # Here we separate components that have already been parameterized and will be restored below.
            clipped = [
                l for l in lines if not (l.startswith("HETATM") and l[17:27] in keep)
            ]
            (d / (pdb + ".pdb")).write_text("\n".join(clipped) + "\n")
            shutil.copyfile(src / pdb / (code + ".cif"), d / (code + ".cif"))
            row = prepare_pose(w, pdb, code, target)
            for name in ("partition_rigid.pdbqt", "receptor.pdbqt"):
                with (d / name).open("a") as f:
                    f.write("\n".join(addons) + "\n")
            row.update(solutes=detail, deposited_pdb_sha256=source["pdb_sha256"])
        except (
            ValueError,
            RuntimeError,
            KeyError,
            subprocess.CalledProcessError,
        ) as exc:
            row.update(status="excluded", reason=str(exc).replace(str(w), "WORK"))
        results.append(row)
        print(row, flush=True)
    (w / "manifest.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
