"""Rebuild external pose RMSDs from deposits, and identify the exclusion contacts."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'cubic_grid/human_panel'))
from verify import reference


def coords(row):
    return [float(row[i:i+8]) for i in (30, 38, 46)]


def atom_id(row):
    return dict(serial=int(row[6:11]), atom=row[12:16].strip(), altloc=row[16].strip(),
                residue=row[17:20].strip(), chain=row[21].strip(), seq=row[22:27].strip(),
                element=row[76:78].strip())


def closest_omitted(path, ligand):
    records = [l for l in path.read_text().splitlines() if l.startswith(('ATOM  ', 'HETATM'))
               and l[76:78].strip() not in ('H', 'D')]
    ligand_rows = [l for l in records if l.startswith('HETATM') and l[17:20].strip() == ligand]
    instance = min(l[21:27] for l in ligand_rows)
    ligand_rows = [l for l in ligand_rows if l[21:27] == instance]
    labels = sorted(set(l[16] for l in ligand_rows) - {' '})
    group = [l for l in ligand_rows if l[16] in {' ', *labels[:1]}]
    query = cKDTree([coords(l) for l in group])
    protein = [l for l in records if l.startswith('ATOM  ')]
    nearest, _ = query.query([coords(l) for l in protein])
    chain = protein[int(np.argmin(nearest))][21]
    candidates = [l for l in records if not (l.startswith('ATOM  ') and l[21] == chain)
                  and l not in ligand_rows and l[17:20].strip() not in ('HOH', 'DOD', 'WAT')]
    distances, index = query.query([coords(l) for l in candidates])
    i = int(np.argmin(distances))
    return dict(receiver_chain=chain, ligand_instance=instance.strip(), distance_A=float(distances[i]),
                ligand_atom=atom_id(group[index[i]]), omitted_atom=atom_id(candidates[i]), altlocs=labels)


def independent_pose(block, mol, exact_reference):
    rows = block.splitlines()
    smiles = next(l.removeprefix('REMARK SMILES ') for l in rows if l.startswith('REMARK SMILES ')
                  and not l.startswith('REMARK SMILES IDX'))
    query = Chem.MolFromSmiles(smiles)
    if query.GetNumAtoms() != mol.GetNumAtoms():
        raise ValueError('Ligand heavy-atom graph differs in size from CCD')
    pairs = []
    for line in rows:
        if line.startswith('REMARK SMILES IDX'):
            pairs.extend(map(int, line.split()[3:]))
    by_id = dict(zip(pairs[::2], pairs[1::2], strict=True))
    xyz = {int(l[6:11]): coords(l) for l in rows if l.startswith(('ATOM  ', 'HETATM'))}
    xyz = np.array([xyz[by_id[i + 1]] for i in range(query.GetNumAtoms())])
    maps = mol.GetSubstructMatches(query, useChirality=True, uniquify=False, maxMatches=10001)
    if not 0 < len(maps) < 10001:
        raise ValueError('Incomplete symmetry enumeration')
    rmsd = min(math.sqrt(np.mean(np.sum((exact_reference[list(m)] - xyz)**2, axis=1))) for m in maps)
    score = float(next(l for l in rows if l.startswith('REMARK VINA RESULT')).split()[3])
    return rmsd, score, len(maps)


def verify(fixtures, runs):
    receipt = {}
    for p, ligand in [('4OCZ', '2RU'), ('6AUM', 'BXV'), ('4R6E', '3JD')]:
        path = fixtures / p
        entry = closest_omitted(path / f'{p}.pdb', ligand)
        state = json.loads((path / 'source.json').read_text())
        assert entry['receiver_chain'] == state['chain']
        assert math.isclose(entry['distance_A'], state['nearest_omitted_nonwater_A'], abs_tol=1e-8)
        assert (entry['distance_A'] >= 8) == (p == '4OCZ')
        receipt[p] = entry
    # The historical reference() resolves the blank / A conformer. Check it covers
    # this actual ligand instance before using it.
    assert receipt['4OCZ']['altlocs'] in ([], ['A'])
    with tempfile.TemporaryDirectory() as temp:
        work = Path(temp)
        (work / 'raw').symlink_to((fixtures / '4OCZ').resolve(), target_is_directory=True)
        _, _, mol = reference(work, '4OCZ', '2RU', receipt['4OCZ']['receiver_chain'])
    ref = mol.GetConformer().GetPositions()
    cells, errors, map_counts = [], [], set()
    for arm in ('official', 'parent', 'stable', 'candidate'):
        folder = runs / arm
        reported = {(r['key'], int(r['rank'])): (float(r['rmsd']), float(r['score']))
                    for r in csv.DictReader((folder / 'ranks.csv').open())}
        cohort = list(csv.DictReader((folder / 'cells.csv').open()))
        assert {r['seed'] for r in cohort} == {'85231', '85237'} and len(cohort) == 2
        for cell in cohort:
            path = folder / 'raw' / (cell['key'] + '.pdbqt')
            assert cell['status'] == 'ok' and hashlib.sha256(path.read_bytes()).hexdigest() == cell['pool_sha256']
            modes = path.read_text().split('ENDMDL')
            assert not modes.pop().strip() and len(modes) == int(cell['modes'])
            rmsds = []
            for rank, block in enumerate(modes, 1):
                rmsd, score, n = independent_pose(block, mol, ref)
                registered, expected_score = reported.pop((cell['key'], rank))
                errors.append(abs(rmsd - registered))
                assert score == expected_score
                rmsds.append(rmsd)
                map_counts.add(n)
            extremes = (rmsds[0], min(rmsds[:5]), min(rmsds), rmsds[0] - min(rmsds))
            for k, v in zip(('first', 'best5', 'best', 'regret'), extremes, strict=True):
                assert abs(float(cell[k]) - v) < 1e-8
            cells.append(dict(arm=arm, seed=cell['seed'], modes=len(rmsds), top1=extremes[0],
                              best5=extremes[1], oracle=extremes[2], regret=extremes[3], sha256=cell['pool_sha256']))
        assert not reported
    assert max(errors) < 1e-8
    return dict(omission=receipt, verified_ranks=len(errors), maximum_rmsd_error_A=max(errors),
                ligand_atoms=mol.GetNumAtoms(), graph_mappings=sorted(map_counts), cells=cells)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixtures', type=Path, help='extracted inputs.tar.xz')
    parser.add_argument('runs', type=Path, help='parent directory of the four arm outputs')
    args = parser.parse_args()
    print(json.dumps(verify(args.fixtures, args.runs), indent=2))
