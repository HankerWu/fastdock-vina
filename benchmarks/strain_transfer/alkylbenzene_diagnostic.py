"""Post-hoc force-field term checks for the two high-affinity T4L alkylbenzenes."""

import argparse
from pathlib import Path

from evaluate import fixed_energy, physical, write
from rdkit import Chem
from rdkit.Chem import AllChem

TERMS = ("Bond", "Angle", "StretchBend", "Torsion", "Oop", "VdW", "Ele")


def term_values(mol, conf_id):
    energies = {}
    for term in TERMS:
        props = physical.mmff_properties(mol)
        for name in TERMS:
            getattr(props, "SetMMFF" + name + "Term")(name == term)
        field = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
        energies[term] = field.CalcEnergy()
    return energies


def analyze(evidence, ligand):
    t4l = evidence / "t4l"
    original = (t4l / "prepared" / f"{ligand}_ligand.pdbqt").read_text()
    proto = physical.Pose(0, 0, 0, original, None, None)
    base, identity = physical.decode_ligand(proto)
    ensemble = Chem.Mol(base)
    ensemble.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = physical.GAS_EMBED_SEED
    params.numThreads = 1
    energies = {}
    for conf in AllChem.EmbedMultipleConfs(
        ensemble, numConfs=physical.GAS_CONFORMERS, params=params
    ):
        energies[conf] = physical.minimized_energy(
            ensemble, physical.mmff_properties(ensemble), conf, ()
        )
    winner = min(energies, key=energies.get)
    reference = term_values(ensemble, winner)
    poses = [("prepared", base)]
    for receptor in ("1L90", "186L"):
        for seed in (271828, 314159):
            path = t4l / "mixed" / f"{receptor}_{ligand}_{seed}_native.pdbqt"
            mol, _ = physical.decode_ligand(physical.read_poses(path)[0])
            poses.append((f"{receptor}_{seed}_top1", mol))
    raw = t4l / "raw" / f"{ligand}-instance.sdf"
    mol = Chem.SDMolSupplier(str(raw), removeHs=True)[0]
    if Chem.MolToSmiles(mol) != identity:
        raise ValueError(f"Instance chemistry differs: {ligand}")
    poses.append(("deposited_crystal", Chem.AddHs(mol, addCoords=True)))
    rows = []
    for source, mol in poses:
        energy = fixed_energy(mol, physical.mmff_properties(mol))
        terms = term_values(mol, mol.GetConformer().GetId())
        deltas = {c: terms[c] - reference[c] for c in TERMS}
        total = energy - energies[winner]
        if abs(sum(deltas.values()) - total) > 1e-5:
            raise ValueError("Term decomposition failed to close")
        rows.append(
            dict(ligand=ligand, source=source, fixed_minus_free=total, **deltas)
        )
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    write(
        args.out,
        [row for ligand in ("184L", "186L") for row in analyze(args.evidence, ligand)],
    )
