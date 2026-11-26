"""Paired 1IEP redocking against deposited STI A 201. Requires RDKit and NumPy."""
import argparse
import csv
from hashlib import sha256
from pathlib import Path
import re
import subprocess

import numpy as np
from rdkit import Chem


BASIC = Path(__file__).resolve().parents[1] / 'examples/basic'
SEEDS = (511, 1013, 2029)
FRESH_SEEDS = (31415, 27182, 16180)


def xyz(line):
    return tuple(float(line[start:start + 8]) for start in (30, 38, 46))


def crystal_reference(pdb):
    crystal = [line for line in pdb.read_text().splitlines()
               if line.startswith('HETATM') and line[17:26] == 'STI A 201'
               and line[76:78].strip() != 'H']
    ligand = (BASIC / '1iep_ligand.pdbqt').read_text().splitlines()
    smiles = next(line.split()[2] for line in ligand if line.startswith('REMARK SMILES '))
    mol = Chem.MolFromSmiles(smiles)
    if len(crystal) != 37 or mol.GetNumAtoms() != 37:
        raise ValueError('Expected the 37 heavy atoms of chain A imatinib')
    by_xyz = {xyz(line): line for line in crystal}
    prepared = {int(line[6:11]): line for line in ligand if line.startswith('ATOM')}
    serials = np.zeros(mol.GetNumAtoms(), dtype=int)
    for line in ligand:
        if line.startswith('REMARK SMILES IDX '):
            indices = list(map(int, line.split()[3:]))
            serials[np.array(indices[::2]) - 1] = indices[1::2]
    ref = []
    pdb_serials = []
    for atom, serial in zip(mol.GetAtoms(), serials):
        # Only the frozen preparation may identify the PDB atom by coordinates.
        # No nearest-neighbor remapping of predicted coordinates is allowed.
        record = by_xyz[xyz(prepared[serial])]
        if atom.GetSymbol() != record[76:78].strip():
            raise ValueError('Ligand element differs from deposited atom')
        ref.append(xyz(record))
        pdb_serials.append(int(record[6:11]))
    if len(set(pdb_serials)) != 37:
        raise ValueError('Input maps ambiguously to the crystal ligand')
    frozen_ids = set(pdb_serials)
    deposited_edges = set()
    for line in pdb.read_text().splitlines():
        if line.startswith('CONECT'):
            ids = [int(line[start:start+5]) for start in range(6, len(line.rstrip()), 5)]
            if ids[0] in frozen_ids:
                deposited_edges.update(tuple(sorted((ids[0], other))) for other in ids[1:] if other in frozen_ids)
    graph_edges = {tuple(sorted((pdb_serials[b.GetBeginAtomIdx()], pdb_serials[b.GetEndAtomIdx()])))
                   for b in mol.GetBonds()}
    if graph_edges != deposited_edges:
        raise ValueError('Chemical embedding differs from PDB connectivity')
    maps = mol.GetSubstructMatches(mol, uniquify=False, useChirality=True)
    return np.array(ref)[np.array(maps)], serials


def dock(binary, seed, options, output, log):
    cmd = [binary, '--receptor', BASIC / '1iep_receptor.pdbqt', '--ligand', BASIC / '1iep_ligand.pdbqt',
           '--config', BASIC / '1iep_receptor.box.txt', '--scoring', 'vina', '--cpu', '1', '--seed', str(seed),
           '--exhaustiveness', '4', '--max_evals', '60000', '--num_modes', '9', '--energy_range', '100',
           '--search_counters', '--out', output]
    cmd.extend(options)
    with log.open('w') as stream:
        subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT, check=True)


def metrics(pool, log, ref, serials):
    rmsds = []
    scores = []
    for block in re.findall(r'^MODEL\s+\d+\s*\n(.*?)^ENDMDL', pool.read_text(), re.M | re.S):
        atoms = {int(l[6:11]): xyz(l) for l in block.splitlines() if l.startswith('ATOM')}
        position = np.array([atoms[s] for s in serials])
        rmsds.append(float(np.sqrt(((ref-position)**2).sum(axis=2).mean(axis=1)).min()))
        scores.append(float(re.search(r'^REMARK VINA RESULT:\s+(\S+)', block, re.M).group(1)))
    counts = re.findall(r'^Search trajectory: .*?hunt_evals=(\d+) .*?second_evals=(\d+)', log.read_text(), re.M)
    if len(counts) != 4:
        raise ValueError('Expected counters for all four trajectories')
    explicit = re.search(r'^Search counters: .*?explicit_evals=(\d+)', log.read_text(), re.M)
    return dict(models=len(scores), first_kcal_mol=scores[0], top1_A=round(rmsds[0], 5),
                best_A=round(min(rmsds), 5), best_rank=1+int(np.argmin(rmsds)),
                grid_evals=sum(int(a)+int(b) for a, b in counts), explicit_evals=int(explicit.group(1)),
                pdbqt_sha256=sha256(pool.read_bytes()).hexdigest())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--mixed', action='store_true', help='Compare native Late vs mixed search, including three new seeds')
    p.add_argument('--old-binary', type=Path, help='Pre-intervention binary for byte-for-byte default checks')
    p.add_argument('--pdb', type=Path, required=True, help='https://files.rcsb.org/download/1IEP.pdb')
    p.add_argument('--out', type=Path, required=True, help='New directory for pools, logs and CSV')
    args = p.parse_args()
    refs, serials = crystal_reference(args.pdb)
    args.out.mkdir(parents=True, exist_ok=False)
    rows = []
    seeds = SEEDS + FRESH_SEEDS if args.mixed else SEEDS
    control = ['--pose_selection', 'late'] if args.mixed else []
    arms = [('native_late', control), ('mixed', ['--mixed_cutoff_search'])] if args.mixed else [
        ('native', []), ('smooth', ['--smooth_cutoff'])]
    for seed in seeds:
        for arm, options in arms:
            pool = args.out / f'{seed}.{arm}.pdbqt'
            log = pool.with_suffix('.log')
            dock(args.binary.resolve(), seed, options, pool, log)
            rows.append(dict(seed=seed, arm=arm, **metrics(pool, log, refs, serials)))
        if args.old_binary:
            old = args.out / f'{seed}.old.pdbqt'
            dock(args.old_binary.resolve(), seed, control, old, old.with_suffix('.log'))
            if old.read_bytes() != (args.out / f'{seed}.{arms[0][0]}.pdbqt').read_bytes():
                raise RuntimeError(f'Default docking changed for seed {seed}')
    with (args.out / 'results.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f'{len(serials)} heavy atoms, {len(refs)} symmetry mappings, {len(rows)} runs')


if __name__ == '__main__':
    main()
