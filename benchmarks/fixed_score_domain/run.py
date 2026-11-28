"""Run the frozen fixed-coordinate domain controls. Python standard library, POSIX."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile
import tempfile
import time

ROOT = Path(__file__).resolve().parents[2]
POSE = ROOT / 'benchmarks/geometry_comparison/pose'
SYSTEMS = ('ABL', 'BPTF')
REFS = ('auto', '0', '-3')
ATOM = ('ATOM  ', 'HETATM')
SCORES = re.compile(r'^\s*(Estimated Free Energy of Binding|\([1-4]\) .*?|Ligand - .*?|Ligand|Flex .*?)\s*:\s*([-+\d.eE]+) \(kcal/mol\)', re.M)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def member(archive, name):
    with tarfile.open(archive, 'r:xz') as tar:
        return tar.extractfile(name).read()


def prepare(work, target):
    if target == 'ABL':
        source = ROOT / 'examples/basic'
        receptor = (source / '1iep_receptor.pdbqt').read_bytes()
        pose = (source / '1iep_ligand.pdbqt').read_text()
    else:
        receptor = member(POSE / 'inputs.tar.xz', '7F5D/receptor.pdbqt')
        pool = member(POSE / 'low_outputs.tar.xz', 'poses/7F5D-7F5D-63419-armijo.pdbqt').decode()
        pose = pool.split('MODEL 1\n', 1)[1].split('ENDMDL', 1)[0]
    work.mkdir()
    (work / 'receptor.pdbqt').write_bytes(receptor)
    rows = pose.splitlines()
    atoms = [row for row in rows if row.startswith(ATOM)]
    coords = [[float(row[i:i+8]) for i in (30, 38, 46)] for row in atoms]
    center = [(min(axis) + max(axis)) / 2 for axis in zip(*coords)]
    for dx in (0, 1, 40, 41):
        translated = [row[:30] + f'{float(row[30:38]) + dx:8.3f}' + row[38:]
                      if row.startswith(ATOM) else row for row in rows]
        (work / f'pose{dx}.pdbqt').write_text('\n'.join(translated) + '\n')
    # These two input ligands contain only ordinary XS heavy types and HD.
    bounded = [float(row[30:38]) for row in atoms if row.split()[-1] != 'HD']
    assert all(row.split()[-1] in {'C', 'A', 'N', 'NA', 'OA', 'SA', 'S', 'HD'} for row in atoms)
    assert min(bounded) + 40 > center[0] + 15
    assert all(max(axis) - min(axis) < 28 for axis in zip(*coords))
    (work / 'malformed.pdbqt').write_text('INVALID\n')
    return center, len(bounded)


def score_blocks(output):
    blocks = []
    for label, value in SCORES.findall(output):
        if label == 'Estimated Free Energy of Binding':
            blocks.append({})
        blocks[-1][label.strip()] = float(value)
    return blocks


def run_cell(binary, work, log, offsets, center, ref, batch):
    args = ['--receptor', 'receptor.pdbqt', '--score_only', '--seed', '20261128', '--cpu', '1']
    poses = [f'pose{n}.pdbqt' if n != 'bad' else 'malformed.pdbqt' for n in offsets]
    if batch:
        args += ['--dir', '.', '--batch', *poses]
    else:
        args += ['--ligand', poses[0]]
    for axis, value in zip('xyz', center):
        args += ['--center_' + axis, f'{value:.6f}', '--size_' + axis, '30']
    if ref != 'auto':
        args += ['--unbound_energy', ref]
    start = time.perf_counter()
    with log.open('wb') as stream:
        proc = subprocess.Popen([binary, *args], cwd=work, stdout=stream, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status)
    elapsed = time.perf_counter() - start
    out = log.read_text()
    return dict(argv=args, exit=proc.returncode, cpu_s=usage.ru_utime + usage.ru_stime,
                wall_s=elapsed, peak_rss_kib=usage.ru_maxrss, blocks=score_blocks(out),
                reported_inputs=re.findall(r'^Scoring ligand: (.*)$', out, re.M),
                outside='outside the grid box' in out, log=log.name)


def evaluate(arms, output):
    receipts = {}
    manifest = {}
    with tempfile.TemporaryDirectory() as raw:
        work = Path(raw)
        for target in SYSTEMS:
            inputs = work / target
            prep_cpu, prep_wall = time.process_time(), time.perf_counter()
            center, bounded = prepare(inputs, target)
            preparation = dict(cpu_s=time.process_time() - prep_cpu, wall_s=time.perf_counter() - prep_wall)
            manifest[target] = dict(center=center, bounded_atoms=bounded, preparation=preparation,
                                    sha256={p.name: digest(p.read_bytes()) for p in sorted(inputs.glob('*.pdbqt'))})
            for arm, binary in arms.items():
                def record(ref, offsets, batch=False, box_delta=0):
                    key = f'{arm}_{target}_{ref}_' + '-'.join(map(str, offsets)) + f'_b{int(batch)}_x{box_delta}'
                    box = [center[0] + box_delta, *center[1:]]
                    receipts[key] = run_cell(binary, inputs, output / f'{key}.log', offsets, box, ref, batch)
                    return receipts[key]
                for ref in REFS:
                    for dx in (0, 1, 40, 41):
                        record(ref, [dx])
                    record(ref, [0], box_delta=1)
                    record(ref, [0, 1], batch=True)
                    if arm == 'candidate':
                        record(ref, [1, 0], batch=True)
                        record(ref, [0, 40, 1], batch=True)
                if arm == 'candidate':
                    record('auto', [0, 'bad', 1], batch=True)
                    for ref in ('nan', 'inf', '-inf'):
                        record(ref, [0])
    return manifest, receipts


def audit(manifest, receipts):
    rows = []
    for target in SYSTEMS:
        def get(arm, ref, poses, batch=0, shift=0):
            return receipts[f'{arm}_{target}_{ref}_{poses}_b{batch}_x{shift}']
        for ref in REFS:
            before = get('stable', ref, 0)
            after = get('candidate', ref, 0)
            assert before['exit'] == after['exit'] == 0
            assert len(before['blocks']) == 1 and len(before['blocks'][0]) == 10
            assert before['blocks'] == after['blocks']
            for arm in ('stable', 'candidate'):
                assert get(arm, ref, 0, shift=1)['blocks'] == before['blocks']
            for dx in (40, 41):
                for arm in ('stable', 'candidate'):
                    result = get(arm, ref, dx)
                    reject = arm == 'candidate' or ref == 'auto'
                    assert (result['exit'] != 0) == reject and result['outside'] == reject
                    assert len(result['blocks']) == int(not reject)
            singles = after['blocks'] + get('candidate', ref, 1)['blocks']
            assert get('stable', ref, '0-1', 1)['blocks'] == []
            for sequence, expected in [('0-1', singles), ('1-0', singles[::-1]), ('0-40-1', singles)]:
                batch = get('candidate', ref, sequence, 1)
                assert batch['blocks'] == expected
                assert batch['reported_inputs'] == [f'pose{d}.pdbqt' for d in sequence.split('-')]
                assert (batch['exit'] != 0) == ('40' in sequence)
            delta = None
            if ref != 'auto':
                receptor = 'Ligand - Receptor'
                delta = (get('stable', ref, 41)['blocks'][0][receptor] -
                         get('stable', ref, 40)['blocks'][0][receptor])
                assert abs(delta - manifest[target]['bounded_atoms'] * 1e6) < 0.01
            rows.append(dict(target=target, reference=ref, bounded_atoms=manifest[target]['bounded_atoms'],
                             interior_score=after['blocks'][0]['Estimated Free Energy of Binding'],
                             old_wall_delta=delta))
        invalid = get('candidate', 'auto', '0-bad-1', 1)
        assert invalid['exit'] and invalid['blocks'] == (get('candidate', 'auto', 0)['blocks'] +
                                                         get('candidate', 'auto', 1)['blocks'])
        for ref in ('nan', 'inf', '-inf'):
            invalid = get('candidate', ref, 0)
            assert invalid['exit'] and invalid['blocks'] == []
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stable', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True)
    arms = {key: getattr(args, key).resolve() for key in ('stable', 'candidate')}
    manifest, receipts = evaluate(arms, args.output)
    meta = dict(inputs=manifest, binaries={key: digest(path.read_bytes()) for key, path in arms.items()})
    (args.output / 'manifest.json').write_text(json.dumps(meta, indent=2) + '\n')
    (args.output / 'runs.json').write_text(json.dumps(receipts, indent=2) + '\n')
    rows = audit(manifest, receipts)
    with (args.output / 'summary.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    print(f'Verified {len(receipts)} real CLI invocations.')


if __name__ == '__main__':
    main()
