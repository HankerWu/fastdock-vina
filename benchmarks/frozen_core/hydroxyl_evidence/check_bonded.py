"""Verify all three analytic torsion profiles against RDKit's independent force field."""
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from bonded import ANALOG, proper_energy, substituents, unit


def check(key):
    smiles, ci, oi = ANALOG[key]
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=43103) == 0
    assert AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s") == 0
    conf = mol.GetConformer()
    xyz = conf.GetPositions()
    c, o = xyz[ci], xyz[oi]
    hyd = next(a.GetIdx() for a in mol.GetAtomWithIdx(oi).GetNeighbors() if a.GetAtomicNum() == 1)
    subs = [a for a in mol.GetAtomWithIdx(ci).GetNeighbors() if a.GetIdx() != oi]
    carbons = [xyz[a.GetIdx()] for a in subs if a.GetSymbol() == "C"]
    mapping = (("CE1","CE2") if key == "TYR" else ("CA","CG2") if key == "THR" else ("CA",))
    coords = dict(zip(mapping, carbons, strict=True))
    coords["CZ" if key == "TYR" else "CB"] = c
    coords["OG1" if key == "THR" else "OG" if key == "SER" else "OH"] = o
    for a, (symbol, vector) in zip(subs, substituents(key, coords, o, c), strict=True):
        assert a.GetSymbol() == symbol
        if symbol == "H":
            xyz[a.GetIdx()] = c + 1.09 * vector
    axis = unit(o-c)
    ang = np.arange(72)*np.pi/36
    vec = xyz[hyd]-o
    proj = np.dot(vec,axis)*axis
    u = vec-proj
    nodes = o + proj + np.outer(np.cos(ang),u) + np.outer(np.sin(ang),np.cross(axis,u))
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    for term in ("Bond", "Angle", "StretchBend", "Oop", "VdW", "Ele"):
        getattr(props, "SetMMFF"+term+"Term")(False)
    field = AllChem.MMFFGetMoleculeForceField(mol, props)
    observed = []
    for node in nodes:
        xyz[hyd] = node
        observed.append(field.CalcEnergy(xyz.flatten().tolist()))
    expected = np.array(observed) - min(observed)
    delta = float(np.max(np.abs(proper_energy(key, coords, nodes) - expected)))
    assert delta < 1.e-8, (key, delta)
    print(f"{key} 72 nodes: max error {delta:.3g} kcal/mol")


if __name__ == "__main__":
    for name in ANALOG:
        check(name)
