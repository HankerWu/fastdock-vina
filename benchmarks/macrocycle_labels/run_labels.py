"""Test physical invariance under renumbering one PDBQT ring-closure identifier."""
import argparse
import csv
import hashlib
import re
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SEED = 458921
CAP = 12000
TRAJECTORY = re.compile(r'hunt_calls=(\d+) hunt_evals=(\d+) second_calls=(\d+) second_evals=(\d+)')
REFINE = re.compile(r' explicit_calls=(\d+) explicit_evals=(\d+)')

def relabel(path, output, label):
    lines = []
    for line in path.read_text().splitlines(keepends=True):
        if line.startswith(('ATOM  ', 'HETATM')):
            line = re.sub(r'(?<=G)0(?=\s*$)', str(label), line)
        lines.append(line)
    output.write_text(''.join(lines))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream-data', required=True, type=Path,
                        help='upstream example/docking_with_macrocycles directory')
    for name in ('official', 'parent', 'fixed', 'full'):
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    a = parser.parse_args()
    data = a.upstream_data / 'solution'
    raw = a.out / 'raw'
    raw.mkdir(parents=True, exist_ok=True)
    receptor = data / 'BACE_1_receptor.pdbqt'
    box = data / 'BACE_1_receptor_vina_box.txt'
    rows = []
    for test in ['single', 'mixed']:
        for label in (range(4) if test == 'single' else (0, 1)):
            input_ligand = raw / f'{test}-{label}.pdbqt'
            original = (data / 'BACE_1_ligand.pdbqt' if test == 'single'
                        else ROOT / 'benchmarks/fixtures/cyclododecane.pdbqt')
            relabel(original, input_ligand, label)
            ligands = ([str(data / 'BACE_1_ligand.pdbqt')] if test == 'mixed' else []) + [str(input_ligand)]
            for name in ('parent', 'fixed', 'full', 'official'):
                if test == 'mixed' and name == 'official':
                    continue
                prefix = f'{test}-{name}-{label}'
                output = raw / (prefix + '.pdbqt')
                base = [str(getattr(a,name).resolve()), '--receptor', str(receptor), '--ligand',
                        *ligands, '--config', str(box), '--cpu', '1', '--seed', str(SEED)]
                if test == 'single':
                    score = subprocess.check_output([*base, '--score_only'], text=True, stderr=subprocess.STDOUT)
                    (raw / (prefix + '-score.log')).write_text(score)
                    fixed_score = float(re.search(r'Estimated Free Energy of Binding\s*:\s*([\-\d.]+)', score)[1])
                else:
                    fixed_score = None
                args = [*base, '--exhaustiveness', '1', '--max_evals', str(CAP),
                        '--num_modes', '5', '--energy_range', '100', '--out', str(output)]
                if name != 'official':
                    args += ['--search_counters']
                start = time.monotonic()
                log = subprocess.check_output(args, text=True, stderr=subprocess.STDOUT)
                elapsed = time.monotonic() - start
                (raw / (prefix + '.log')).write_text(log)
                coords = '\n'.join(line[:54] for line in output.read_text().splitlines()
                                   if line.startswith(('ATOM  ', 'HETATM')))
                counts = TRAJECTORY.search(log)
                refine = REFINE.search(log)
                if name != "official" and (not counts or not refine):
                    raise RuntimeError(f"Missing evaluation counters in {prefix}")
                rows.append(dict(test=test, engine=name, label=label, fixed_score=fixed_score,
                                 coordinate_hash=hashlib.sha256(coords.encode()).hexdigest(),
                                 hunt_calls=counts[1] if counts else None, hunt_evals=counts[2] if counts else None,
                                 second_calls=counts[3] if counts else None, second_evals=counts[4] if counts else None,
                                 explicit_calls=refine[1] if refine else None, explicit_evals=refine[2] if refine else None,
                                 wall_s=round(elapsed,4)))
    with (a.out / 'runs.csv').open('w') as handle:
        writer = csv.DictWriter(handle, rows[0])
        writer.writeheader()
        writer.writerows(rows)
    for test in ('single','mixed'):
        fixed = [r for r in rows if r['engine'] == 'fixed' and r['test'] == test]
        for key in ('coordinate_hash', 'fixed_score', 'hunt_evals', 'second_evals', 'explicit_evals'):
            if len({r[key] for r in fixed}) != 1:
                raise RuntimeError(f'{test} changed {key} under renaming')

if __name__ == '__main__':
    main()
