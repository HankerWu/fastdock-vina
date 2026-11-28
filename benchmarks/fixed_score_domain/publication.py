"""Replay native and late published pools on their original ligand graph."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FLAGGED = {(85231, r) for r in (5, 7, 8)} | {(85237, 8)}


def run(cmd, log):
    start = time.perf_counter()
    with log.open('w') as stream:
        p = subprocess.Popen(cmd, stdout=stream, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(p.pid, 0)
        p.returncode = os.waitstatus_to_exitcode(status)
    return dict(log=log.name, argv=list(map(str, cmd)), exit=p.returncode, cpu_s=usage.ru_stime+usage.ru_utime,
                wall_s=time.perf_counter()-start, peak_rss_kib=usage.ru_maxrss)


def energy(text, label):
    line = next(x for x in text.splitlines() if x.startswith(label))
    return float(line.split(':', 1)[1].split()[0])


def coords(row):
    return [float(row[i:i+8]) for i in (30, 38, 46)]


def csvsave(path, rows):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('before', 'after', 'output'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    work = a.output.resolve()
    work.mkdir(parents=True)
    raw = work / 'raw'
    raw.mkdir()
    with tarfile.open(ROOT / 'benchmarks/cumulative/external/inputs.tar.xz') as tar:
        tar.extractall(work / 'fixtures', filter='data')
    case = work / 'fixtures/4OCZ'
    meta = json.loads((case / 'metadata.json').read_text())
    center = np.array(meta['center'])
    reference = np.array(meta['mapped'])
    base = ROOT / 'examples/basic/1iep'
    costs, cells, poses, points, logged_counts = [], [], [], {}, {}
    binary_hash = {}
    for arm, binary in [('before', a.before.resolve()), ('after', a.after.resolve())]:
        binary_hash[arm] = hashlib.sha256(binary.read_bytes()).hexdigest()
        for config, seeds in [('native', [85231,85237]), ('late', [85231,85237]), ('abl', [85231])]:
            for seed in seeds:
                key = f'{arm}-{config}-{seed}'
                rec, lig = case / 'receptor.pdbqt', case / 'ligand.pdbqt'
                box = [item for axis, value in zip('xyz', center) for item in
                       ('--size_' + axis, '24', '--center_' + axis, str(value))]
                if config == 'abl':
                    rec, lig = [Path(str(base) + suffix) for suffix in ('_receptor.pdbqt', '_ligand.pdbqt')]
                    box = ['--config', Path(str(base) + '_receptor.box.txt')]
                cmd = [binary, '--ligand', lig, '--receptor', rec, '--cpu','1','--seed',str(seed),
                       '--spacing','0.375', *box]
                out = raw / (key + '.pdbqt')
                log = raw / (key + '.log')
                flags = ['--pose_selection','late'] if config == 'late' else []
                c = run([*cmd, '--exhaustiveness','3','--max_evals','12000','--num_modes','9',
                         '--energy_range','100','--search_counters',*flags,'--out',out], log)
                costs.append(c)
                if c['exit']:
                    raise ValueError((key,log.read_text()))
                counts = [l for l in log.read_text().splitlines()
                          if l.startswith(('Search trajectory:', 'Search counters:'))]
                assert len(counts) == 4
                logged_counts[arm, config, seed] = counts
                cell = dict(arm=arm, config=config, seed=seed, hash=hashlib.sha256(out.read_bytes()).hexdigest())
                blocks = out.read_text().split('ENDMDL')[:-1]
                cell['modes'] = len(blocks)
                cells.append(cell)
                if config != 'native':
                    continue
                for rank, block in enumerate(blocks, 1):
                    pose = raw / f'{key}-{rank}.pdbqt'
                    rows = [r for r in block.strip().splitlines() if not r.startswith('MODEL')]
                    pose.write_text('\n'.join(rows) + '\n')
                    log = pose.with_suffix('.log')
                    c = run([*cmd, '--score_only','--pose_coordinates',pose,'--unbound_energy','0'],log)
                    costs.append(c)
                    atoms = {int(l[6:11]):coords(l) for l in rows if l.startswith(('ATOM  ', 'HETATM'))}
                    heavy = np.asarray([atoms[i] for i in meta['serial']])
                    points[arm, seed, rank] = heavy
                    rmsd = np.sqrt(np.mean(np.sum((reference - heavy)**2, axis=2), axis=1)).min()
                    exceed = np.maximum(0, np.abs(heavy - center) - 12)
                    archived = energy(block, 'REMARK INTER + INTRA:')
                    score = log.read_text()
                    if c['exit'] == 0:
                        rescored = energy(score, '(1) Final Intermolecular') + energy(score, '(2) Final Total Internal')
                        delta = rescored - archived
                    else:
                        rescored = delta = None
                    poses.append(dict(arm=arm, seed=seed, rank=rank, flagged=(seed, rank) in FLAGGED,
                                      scored=c['exit']==0, archived=archived, fixed_raw=rescored,
                                      residual=delta, excess_A=exceed.max(), wall_sum=1e6*exceed.sum(),
                                      crystal_rmsd=rmsd, shift_rmsd=0., max_atom_shift=0.))
    for row in poses:
        prior = points['before', row['seed'], row['rank']]
        current = points[row['arm'], row['seed'], row['rank']]
        distances = np.linalg.norm(current-prior, axis=1)
        row.update(shift_rmsd=np.sqrt(np.mean(distances**2)), max_atom_shift=distances.max())
    for cell in cells:
        peer = next(r for r in cells if r['arm']=='before' and r['config']==cell['config'] and r['seed']==cell['seed'])
        cell['same_pool_bytes'] = cell['hash'] == peer['hash']
        cell['same_gradient_counts'] = (logged_counts[cell['arm'], cell['config'], cell['seed']]
                                         == logged_counts['before', cell['config'], cell['seed']])
    csvsave(work / 'modes.csv', poses)
    csvsave(work / 'runs.csv', cells)
    (work / 'costs.json').write_text(json.dumps(costs,indent=2)+'\n')
    (work / 'binaries.json').write_text(json.dumps(binary_hash,indent=2)+'\n')
    assert len(poses) == 36 and all(p['scored'] for p in poses)
    assert all(c['same_gradient_counts'] for c in cells)
    assert all(c['same_pool_bytes'] for c in cells if c['config']=='late')
    assert all(p['wall_sum'] < 1 and abs(p['residual']) < 0.5 for p in poses if p['arm']=='after' and p['flagged'])


if __name__ == '__main__':
    main()
