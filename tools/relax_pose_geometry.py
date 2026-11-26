"""Relax internal ligand geometry with harmonic tethers to a Vina pose in its receptor frame."""

import argparse
import math
from pathlib import Path

from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import AllChem


def first_pose(pool):
    mols = RDKitMolCreate.from_pdbqt_mol(
        PDBQTMolecule.from_file(str(pool), skip_typing=True)
    )
    if len(mols) != 1 or mols[0] is None:
        raise ValueError("Expected one mapped Meeko ligand in the native pose pool")
    mol = mols[0]
    top = Chem.Mol(mol)
    top.RemoveAllConformers()
    top.AddConformer(Chem.Conformer(mol.GetConformer()), assignId=True)
    return Chem.AddHs(top, addCoords=True)


def relax(mol, stiffness=10.0):
    if stiffness < 0 or not math.isfinite(stiffness):
        raise ValueError("stiffness must be finite and nonnegative")
    if mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D():
        raise ValueError("Expected one 3D ligand pose")
    moved = Chem.Mol(mol)
    params = AllChem.MMFFGetMoleculeProperties(moved, mmffVariant="MMFF94s")
    if params is None:
        raise ValueError("MMFF94s parameters unavailable for this ligand")
    ff = AllChem.MMFFGetMoleculeForceField(moved, params, nonBondedThresh=100.0)
    for atom in moved.GetAtoms():
        if atom.GetAtomicNum() > 1 and stiffness:
            ff.MMFFAddPositionConstraint(atom.GetIdx(), 0.0, stiffness)
    ff.Initialize()
    initial = ff.CalcEnergy()
    status = ff.Minimize(maxIts=2000, forceTol=1e-4, energyTol=1e-7)
    final = ff.CalcEnergy()
    if status or not math.isfinite(final) or final > initial + 1e-6:
        raise ValueError("Tethered ligand minimization did not converge")
    return moved


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("native_pool", type=Path)
    p.add_argument("output_sdf", type=Path)
    p.add_argument(
        "--stiffness", type=float, default=10.0, help="kcal/mol/Å² for each heavy atom"
    )
    a = p.parse_args()
    mol = relax(first_pose(a.native_pool), a.stiffness)
    mol.SetProp("pose_source", "native rank 1")
    mol.SetDoubleProp("heavy_atom_tether_kcal_mol_A2", a.stiffness)
    # Scores stored on a docking input describe the original coordinates.
    # The exported record carries provenance and the geometry operation only.
    for prop in list(mol.GetPropNames()):
        if prop not in ("pose_source", "heavy_atom_tether_kcal_mol_A2"):
            mol.ClearProp(prop)
    with Chem.SDWriter(str(a.output_sdf)) as w:
        w.write(mol)


if __name__ == "__main__":
    main()
