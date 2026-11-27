"""Add the proper torsion about the alcohol or phenol C–O bond to a dry receptor prior."""
from collections import defaultdict
from functools import lru_cache

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Atom order is explicit for the tiny reference fragments. RDKit supplies all
# X–C–O–H parameters, including X=H; MMFF94s uses kcal/mol.
ANALOG = {"SER": ("CCO", 1, 2), "THR": ("CC(O)C", 1, 2), "TYR": ("Oc1ccccc1", 1, 0)}


def unit(x):
    length = np.linalg.norm(x)
    if length < 1.e-8:
        raise ValueError("Undefined local hydroxyl frame")
    return x / length


@lru_cache(maxsize=None)
def torsions(residue):
    smiles, c, o = ANALOG[residue]
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    hydro = next(a.GetIdx() for a in mol.GetAtomWithIdx(o).GetNeighbors() if a.GetAtomicNum() == 1)
    terms = []
    for atom in mol.GetAtomWithIdx(c).GetNeighbors():
        if atom.GetIdx() == o:
            continue
        params = props.GetMMFFTorsionParams(mol, atom.GetIdx(), c, o, hydro)
        if params is None:
            raise ValueError("Missing proper torsion parameters")
        terms.append((atom.GetSymbol(), params[1:]))
    return terms


def substituents(residue, coords, omega, pivot):
    if residue == "TYR":
        return [("C", coords[k] - pivot) for k in ("CE1", "CE2")]
    to_alpha = unit(coords["CA"] - pivot)
    to_oxygen = unit(omega - pivot)
    if residue == "THR":
        methyl = unit(coords["CG2"] - pivot)
        return [("C", to_alpha), ("C", methyl), ("H", -unit(to_alpha + methyl + to_oxygen))]
    bisector = unit(to_alpha + to_oxygen)
    cross = unit(np.cross(to_alpha, to_oxygen))
    # Choose symmetric H directions with tetrahedral cosines to both heavy neighbors.
    axial = -1 / (3 * np.dot(bisector, to_alpha))
    if abs(axial) >= 1:
        raise ValueError("Incompatible tetrahedral geometry")
    offset = np.sqrt(1 - axial**2) * cross
    return [("C", to_alpha), ("H", axial * bisector + offset), ("H", axial * bisector - offset)]


def proper_energy(residue, coords, points):
    pivot, oxygen = ((coords["CZ"], coords["OH"]) if residue == "TYR" else
                     (coords["CB"], coords["OG1" if residue == "THR" else "OG"]))
    axis = unit(oxygen - pivot)
    end = points - oxygen
    end = end - np.outer(end @ axis, axis)
    end /= np.linalg.norm(end, axis=1)[:, None]
    energy = np.zeros(len(points))
    groups = substituents(residue, coords, oxygen, pivot)
    for (symbol, params), (name, vector) in zip(torsions(residue), groups, strict=True):
        if name != symbol:
            raise ValueError("Analogue mapping disagrees with receptor neighbor")
        c = np.clip(end @ unit(vector - axis * np.dot(axis, vector)), -1, 1)
        v1, v2, v3 = params
        energy += .5 * (v1*(1+c) + v2*(2-2*c*c) + v3*(1 + 4*c**3 - 3*c))
    return energy - energy.min()


def add_bonded_prior(receptor):
    """Return a shallow copy; share geometry and replace only the prior energies."""
    from copy import copy
    residues = defaultdict(dict)
    for a in receptor.atoms.values():
        residues[a.residue][a.name] = a.xyz
    result = copy(receptor)
    result.rotors = {}
    for serial, (nodes, energies) in receptor.rotors.items():
        residue = receptor.atoms[serial].residue
        bonded = proper_energy(residue[:3], residues[residue], nodes)
        result.rotors[serial] = nodes, energies + bonded
    return result
