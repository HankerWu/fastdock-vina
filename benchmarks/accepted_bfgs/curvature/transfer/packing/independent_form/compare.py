"""Evaluate archived VH032 poses against the independent monoclinic 9GIO crystal."""
import argparse
import gzip
import hashlib
import itertools
import json
import sys
import tarfile
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem, RDLogger
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent.parent
sys.path.insert(0, str(STUDY / "crystal_copies"))
from replay import CHAINS, PDBS, read_structure, register, rmsd, transform, write
sys.path.insert(0, str(HERE.parent))
from measure import Crystal, heavy_modes


def crystal_data(path, smiles):
    structure = gemmi.read_structure(str(path))
    structure.setup_entities()
    chain = structure[0]["C"]
    ligand = next(r for r in chain if r.name == "3JF")
    atoms = [a for a in ligand if not a.element.is_hydrogen]
    assert len(atoms) == 33 and all(a.altloc == "\x00" and a.occ == 1 for a in atoms)
    xyz = np.array([list(a.pos) for a in atoms])
    block = gemmi.cif.read(str(path)).sole_block()
    ids = {a.name: i for i, a in enumerate(atoms)}
    mol = Chem.RWMol()
    conformer = Chem.Conformer(len(atoms))
    for i, a in enumerate(atoms):
        mol.AddAtom(Chem.Atom(a.element.name))
        conformer.SetAtomPosition(i, tuple(xyz[i]))
    types = {"sing": Chem.BondType.SINGLE, "doub": Chem.BondType.DOUBLE,
             "trip": Chem.BondType.TRIPLE}
    for comp, a, b, order in block.find("_chem_comp_bond.", [
            "comp_id", "atom_id_1", "atom_id_2", "value_order"]):
        if comp == "3JF" and a in ids and b in ids:
            mol.AddBond(ids[a], ids[b], types[order])
    mol.AddConformer(conformer)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
    key = Chem.MolToInchiKey(mol)
    assert key == "GFVIEZBZIUKYOG-SVFBPWRDSA-N"
    matches = mol.GetSubstructMatches(Chem.MolFromSmiles(smiles), useChirality=True,
                                      uniquify=False, maxMatches=10001)
    assert len(matches) == 12
    ca = {}
    for res in chain:
        if res.entity_type != gemmi.EntityType.Polymer:
            continue
        options = [a for a in res if a.name == "CA" and a.occ > 0]
        if options:
            key_ca = (f"{res.seqid.num:4}{res.seqid.icode}", res.name)
            ca[key_ca] = [(a.altloc, list(a.pos)) for a in options]
    q = next(r for r in chain if r.name == "A1IMD")
    adduct = np.array([list(a.pos) for a in q if a.occ > 0 and not a.element.is_hydrogen])
    info = dict(ccd="3JF", inchikey=key, heavy_atoms=len(atoms), symmetry_maps=len(matches),
                space_group=structure.spacegroup_hm, resolution_A=structure.resolution,
                closest_adduct_atom_A=float(cdist(xyz, adduct).min()),
                cif_sha256=hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest())
    return xyz[list(matches)], ca, info


def compare(work, out):
    RDLogger.DisableLog("rdApp.warning")
    with tarfile.open(STUDY / "evidence.tar.xz") as archive:
        old = {p: read_structure(archive, p) for p in PDBS}
    q, r = PDBS
    meta = json.loads((work / q / "metadata.json").read_text())
    new, ca, info = crystal_data(HERE / "9gio.cif.gz", meta["smiles"])
    old_refs, _ = register(old, q, r, "pocket")
    keys = sorted(old[q][3] & old[r][3])
    assert len(keys) == 55 and all(len(ca[k]) == 1 for k in keys)
    a, b = np.array([ca[k][0][1] for k in keys]), np.array([old[r][2]["C"][k] for k in keys])
    new_ref = transform(a, b, new)
    info.update(ca_count=len(keys), pocket_rmsd_A=rmsd(transform(a, b, a), b))
    copies = {ch: rmsd(old_refs[ch][0], new_ref) for ch in CHAINS}
    source_ca = np.array([old[q][2]["C"][k] for k in keys])
    info["direct_9GIO_to_4W9H_C_A"] = rmsd(transform(a, source_ca, new), old[q][1]["C"][0])
    entries = [("deposited_9GIO", "", "", 1, new_ref[0])]
    for ch in CHAINS:
        entries.append(("deposited_4W9H", ch, "", 1, old_refs[ch][0]))
    for group, files in (("blind", (work/"runs").glob(f"{q}_{r}_*pdbqt")),
                         ("explicit_graph", work.glob(f"{q}_{r}_*pdbqt"))):
        for f in sorted(files):
            _, _, seed, arm = f.stem.split("_", 3)
            for rank, (xyz, _) in enumerate(heavy_modes(f, meta["serial"]), 1):
                entries.append((group, arm, seed, rank, xyz))
    witness, _ = next(heavy_modes(work/"witness"/f"{q}_fixed_local.pdbqt", meta["serial"]))
    entries.append(("witness", "fixed", "", 1, witness))
    in_source = [transform(b, a, e[-1]) for e in entries]
    env = Crystal(HERE / "9gio.cif.gz", np.concatenate(in_source))
    # The biological assembly is the only trimer in this asymmetric unit.
    assert set(env.trees) == {"image"}
    rows = []
    for (group, arm, seed, rank, pose), xyz in zip(entries, in_source):
        distances = env.measure(xyz, meta["serial"])
        rows.append(dict(group=group, arm=arm, seed=seed, rank=rank,
                         rmsd_9GIO_A=rmsd(pose, new_ref), rmsd_4W9H_A=rmsd(pose, old_refs["C"]),
                         closest_9GIO_mate_A=distances["min_A"], pairs_lt_2p3=distances["clashes"],
                         ligand_serial=distances["image_ligand_serial"], mate=distances["image_atom"],
                         mate_op=distances["image_op"]))
    blind = [v for v in rows if v["group"] == "blind"]
    info.update(copies_to_9GIO_A=copies, n_blind=len(blind),
                max_blind_reference_change_A=max(abs(v["rmsd_9GIO_A"] - v["rmsd_4W9H_A"]) for v in blind))
    assert info["max_blind_reference_change_A"] <= copies["C"] + 1e-10
    summary = []
    for (arm, seed), pool in itertools.groupby(sorted(blind, key=lambda v: (v["arm"], v["seed"])),
                                               lambda v: (v["arm"], v["seed"])):
        obs = list(pool)
        assert obs[0]["rank"] == 1
        summary.append(dict(arm=arm, seed=seed, modes=len(obs),
                            top1_4W9H_A=obs[0]["rmsd_4W9H_A"], top1_9GIO_A=obs[0]["rmsd_9GIO_A"],
                            oracle_4W9H_A=min(v["rmsd_4W9H_A"] for v in obs),
                            oracle_9GIO_A=min(v["rmsd_9GIO_A"] for v in obs)))
    out.mkdir(parents=True, exist_ok=True)
    write(out/"forms.csv", rows)
    write(out/"pools.csv", summary)
    (out/"comparison.json").write_text(json.dumps(info, indent=2) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument("out", type=Path)
    args = ap.parse_args()
    compare(args.work, args.out)
