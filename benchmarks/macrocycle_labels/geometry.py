"""Audit E4Y in the 5QCT chain-A reference frame and the output closure geometry."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from rdkit import Chem


def xyz(line):
    return np.array([float(line[j:j+8]) for j in (30,38,46)])


def alpha_carbons(path):
    return {line[17:27]: xyz(line) for line in path.read_text().splitlines()
            if line.startswith('ATOM  ') and line[12:16].strip() == 'CA' and line[21] == 'A'
            and line[16] in (' ', 'A')}


def coord_blocks(path):
    current = {}
    for line in path.read_text().splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            current[int(line[6:11])] = (line.split()[-1], xyz(line))
        elif line.startswith('ENDMDL'):
            yield current
            current = {}
    if current:
        yield current


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream-data', type=Path, required=True)
    parser.add_argument('--structure', type=Path, required=True, help='5QCT.pdb from wwPDB')
    parser.add_argument('--out', type=Path, required=True, help='output directory from run.py')
    a = parser.parse_args()
    solution = a.upstream_data / 'solution'
    query_ca = alpha_carbons(solution/'BACE_1_receptor.pdbqt')
    crystal_ca = alpha_carbons(a.structure)
    names = sorted(query_ca.keys() & crystal_ca.keys())
    query = np.array([query_ca[k] for k in names])
    crystal = np.array([crystal_ca[k] for k in names])
    # One rigid receptor fit. Never align a predicted ligand.
    u, _, vt = np.linalg.svd((crystal-crystal.mean(0)).T @ (query-query.mean(0)))
    u[:, -1] *= np.linalg.det(u @ vt)
    rot = u @ vt
    def align(x):
        return (x - crystal.mean(0)) @ rot + query.mean(0)

    mol2 = a.upstream_data / 'data/BACE_1_ligand.mol2'
    mol = Chem.MolFromMol2File(str(mol2), removeHs=True)
    coords = mol.GetConformer().GetPositions()
    maps = np.array(mol.GetSubstructMatches(mol, uniquify=False,useChirality=True,maxMatches=10000))
    if len(maps) == 10000:
        raise RuntimeError("Truncated graph automorphisms")
    deposited = {line[12:16].strip(): xyz(line) for line in a.structure.read_text().splitlines()
                 if line.startswith('HETATM') and line[17:20]=='E4Y' and line[21:26]=='A 401'}
    ref = align(np.array([deposited[atom.GetProp('_TriposAtomName')] for atom in mol.GetAtoms()]))
    inp = next(coord_blocks(solution / 'BACE_1_ligand.pdbqt'))
    order = []
    for coord in coords:
        matches = [i for i, (typ, pos) in inp.items()
                   if typ not in ('HD', 'G0') and np.linalg.norm(pos - coord) < .002]
        if len(matches) != 1:
            raise ValueError("Expected one original MOL2 heavy atom at the PDBQT coordinate")
        order.append(matches[0])
    if len(set(order)) != len(coords):
        raise ValueError("Incomplete atom correspondence")
    def rmsd(pos):
        return float(np.min(np.sqrt(np.mean(np.sum((pos[maps] - ref)**2,axis=2),axis=1))))
    atom_a, atom_b = [i for i,(typ,pos) in inp.items() if typ=='CG0']
    dummies = [i for i,(typ,pos) in inp.items() if typ=='G0']
    pairs = [(i, min(dummies,key=lambda j:np.linalg.norm(inp[i][1]-inp[j][1]))) for i in (atom_a,atom_b)]
    if any(np.linalg.norm(inp[i][1]-inp[j][1])>0.002 for i,j in pairs):
        raise ValueError("Pseudoatoms must start on their opposite closure carbon")
    info = dict(pdb='5QCT', chain='A', ligand='E4Y', common_receptor_CAs=len(names),
                CA_fit_rmsd=float(np.sqrt(np.mean(np.sum((align(crystal)-query)**2,axis=1)))),
                input_rmsd=rmsd(coords), chemical_heavy_atoms=len(order), graph_automorphisms=len(maps))
    (a.out/'registration.json').write_text(json.dumps(info,indent=2)+'\n')
    results = []
    for posefile in sorted((a.out/'raw').glob('single-*-?.pdbqt')):
        _, engine,label = posefile.stem.split('-')
        for rank,pose in enumerate(coord_blocks(posefile),1):
            pos = np.array([pose[i][1] for i in order])
            results.append(dict(engine=engine,label=label,rank=rank,rmsd=round(rmsd(pos),6),
                                bond_length=round(np.linalg.norm(pose[atom_a][1]-pose[atom_b][1]),6),
                                max_dummy_residual=round(max(np.linalg.norm(pose[i][1]-pose[j][1]) for i,j in pairs),6)))
    with (a.out/'geometry.csv').open('w') as handle:
        writer = csv.DictWriter(handle,results[0])
        writer.writeheader()
        writer.writerows(results)

if __name__ == '__main__':
    main()
