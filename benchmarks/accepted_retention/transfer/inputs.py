"""Prepare independent human pocket pairs with their resolved crystal buffer components."""

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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pair_gradient/pose"))
import paired as reference
from sources import crystal_contacts

base = reference.base
CASES = {
    "4N1U": ("NUDT1", "2GE", "A"),
    "4C9X": ("NUDT1", "VHS", "A"),
    "5NRA": ("CHIT1", "95K", "A"),
    "6ZE8": ("CHIT1", "QGB", "C"),
}
SEEDS = (214087, 214091)
CAPS = (12000, 24000)


def addon(lines, ccd, first_id):
    """Type crystal coordinates by CCD covalent graph. Polar hydrogens are generated."""
    comp = gemmi.cif.read_file(str(ccd)).sole_block()
    byname = {l[12:16].strip(): l for l in lines if l[16] in " A"}
    atoms = comp.get_mmcif_category("_chem_comp_atom.")
    graph = Chem.RWMol()
    indices = {}
    for name, e, q in zip(atoms["atom_id"], atoms["type_symbol"], atoms["charge"]):
        if e in ("H", "D"):
            continue
        a = Chem.Atom(e.title())
        a.SetFormalCharge(int(q))
        indices[name] = graph.AddAtom(a)
    assert indices.keys() == byname.keys()
    bonds = comp.get_mmcif_category("_chem_comp_bond.")
    for a, b, kind in zip(
        bonds.get("atom_id_1", []),
        bonds.get("atom_id_2", []),
        bonds.get("value_order", []),
    ):
        if a in indices and b in indices:
            graph.AddBond(
                indices[a],
                indices[b],
                {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE}[kind],
            )
    mol = graph.GetMol()
    Chem.SanitizeMol(mol)
    conf = Chem.Conformer(len(indices))
    for name, i in indices.items():
        conf.SetAtomPosition(i, base.coords(byname[name]))
    mol.AddConformer(conf)
    if len(indices) == 1:
        # Meeko's ligand writer does not need to construct a torsion tree for a single chloride.
        atom = mol.GetAtomWithIdx(0)
        assert atom.GetSymbol() == "Cl" and atom.GetFormalCharge() == -1
        line = next(iter(byname.values()))
        return [line[:6] + f"{first_id:5d}" + line[11:66] + "    -1.000 Cl"]
    h = Chem.AddHs(mol, addCoords=True)
    (setup,) = MoleculePreparation().prepare(h)
    block, ok, err = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(err)
    out = []
    residue = lines[0][17:27]
    for line in block.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            out.append(
                "HETATM"
                + f"{first_id + len(out):5d}"
                + line[11:17]
                + residue
                + line[27:]
            )
    return out


def prepare(work, pdb):
    target, code, ch = CASES[pdb]
    out = work / pdb
    raw = work / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    path, ccd = raw / (pdb + ".pdb"), raw / (code + ".cif")
    psha = base.fetch(path, f"https://files.rcsb.org/download/{pdb}.pdb")
    csha = base.fetch(ccd, f"https://files.rcsb.org/ligands/download/{code}.cif")
    lines = path.read_text().splitlines()
    atom_lines = [
        s
        for s in lines
        if s.startswith(("ATOM  ", "HETATM")) and s[76:78].strip() not in {"H", "D"}
    ]
    chosen = [s[:21] + "A" + s[22:] for s in atom_lines if s[21] == ch]
    refs = base.dep_records(chosen, code)
    crystal = np.array(list(refs[0].values()))
    q = base.ligand(ccd, refs, out / "ligand.pdbqt")
    ca = {
        s[22:27].strip() + ":" + s[17:20]: base.coords(s)
        for s in chosen
        if s.startswith("ATOM  ") and s[12:16].strip() == "CA" and s[16] in " A"
    }
    # Retain complete groups touching the pocket, fail on unmodeled molecular neighbors.
    groups = {}
    for s in atom_lines:
        if s[17:20] == "HOH" or (
            s[21] == ch and (s.startswith("ATOM") or s[17:20].strip() == code)
        ):
            continue
        groups.setdefault(s[17:27], []).append(s)
    distances = {
        k: float(
            np.min([np.linalg.norm(crystal - base.coords(s), axis=1).min() for s in ss])
        )
        for k, ss in groups.items()
    }
    keep = [k for k, dist in distances.items() if dist < 8]
    if any(k[:3].strip() not in {"SO4", "CL", "GOL"} or k[4] != ch for k in keep):
        raise ValueError((pdb, "unmodeled near-pocket neighbor", keep))
    # Ligand and buffer atom coordinates never enter protein completion.
    rec = base.receptor([s for s in chosen if s.startswith("ATOM")], out, pdb, crystal)
    if (
        rec["closest_repaired_ligand_distance"] is not None
        and rec["closest_repaired_ligand_distance"] < 8
    ):
        raise ValueError((pdb, rec["closest_repaired_ligand_distance"]))
    retained = []
    receptor = out / "receptor.pdbqt"
    with receptor.open("a") as f:
        next_id = 8000
        for k in sorted(keep):
            comp = k[:3].strip()
            filename = raw / (comp + ".cif")
            digest = base.fetch(
                filename, f"https://files.rcsb.org/ligands/download/{comp}.cif"
            )
            new = addon(groups[k], filename, next_id)
            next_id += len(new)
            f.write("\n".join(new) + "\n")
            retained.append(
                {
                    "group": k,
                    "distance_A": distances[k],
                    "ccd_sha256": digest,
                    "atoms_written": len(new),
                    "atom_types": [s.split()[-1] for s in new],
                    "minimum_occupancy": min(float(s[54:60]) for s in groups[k]),
                }
            )
    omitted = sorted((v, k) for k, v in distances.items() if k not in keep)
    structure = gemmi.read_structure(str(path))
    metadata = dict(
        pdb=pdb,
        target=target,
        chain=ch,
        ligand=code,
        pdb_sha256=psha,
        ccd_sha256=csha,
        resolution_A=structure.resolution,
        ca=ca,
        center=crystal.mean(axis=0).tolist(),
        receptor_status="ready",
        retained_additives=retained,
        nearest_omitted=omitted[0] if omitted else None,
        ligand_minimum_occupancy=min(
            float(s[54:60]) for s in chosen if s[17:20].strip() == code
        ),
        receptor_sha256=hashlib.sha256(receptor.read_bytes()).hexdigest(),
        **crystal_contacts(structure, ch, code),
        **q,
        **rec,
    )
    (out / "preparation.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(pdb, "ready", retained, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("pdb", nargs="*")
    args = p.parse_args()
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    for key in args.pdb or CASES:
        prepare(args.work, key)
