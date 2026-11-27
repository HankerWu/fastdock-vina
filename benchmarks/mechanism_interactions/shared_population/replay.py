"""Test simultaneous two-state population constraints against within-assay Ki or Kd."""
import argparse
import csv
import hashlib
import itertools
import json
import lzma
import math
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

from scipy.stats import kendalltau, spearmanr

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
KINASE = BENCH / 'kinase_state_affinity'
sys.path.insert(0, str(KINASE / 'implementation'))
from ensemble_affinity import R, pair_sensitivity

RT = R * 298.15
HALF = Fraction(1, 2)


def read(path):
    opener = lzma.open if path.suffix == '.xz' else open
    with opener(path, 'rt', newline='') as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def collect():
    paths = [KINASE / 'panel.csv', KINASE / 'results/cells.csv',
             KINASE / 'independent-seeds/results/cells.csv',
             BENCH / 'frozen_core/strain_rotor_interaction/first_poses.csv.xz']
    rows = defaultdict(list)
    labels = {(r['target'], r['cid']): r['kd_nm'] for r in read(paths[0])}
    for path in paths[1:3]:
        for r in read(path):
            t, ligand = r['target'], r['ligand']
            for column, ref in [('native', 'pool'), ('energy', 'prepared')]:
                key = (t, 'Davis-' + t, 'Kd', ref)
                rows[key].append((ligand, r['state'], r['replicate'], labels[t, ligand], r[column]))
    for r in read(paths[3]):
        panel = r['panel']
        if panel not in {'native_12k', 'fxa_repaired', 'parp_spr', 'transfer_archived', 'transfer_active_torsions'}:
            continue
        reference = 'active_rotor_control' if panel == 'transfer_active_torsions' else 'pool'
        key = (r['target'], r['assay'], r['endpoint'], reference)
        rows[key].append((r['compound'], r['receptor'], r['seed'], r['measured_nM'], r['native']))
    return {key: records for key, records in rows.items()
            if len({r[1] for r in records}) == 2}, paths


@dataclass(frozen=True)
class Ligand:
    name: str
    nm: Decimal
    energy: tuple
    activity: tuple


def prepare(records):
    by_ligand = defaultdict(dict)
    measurements = {}
    for ligand, state, seed, nm, energy in records:
        nm = Decimal(nm)
        energy = Decimal(energy)
        if not nm.is_finite() or nm <= 0 or not energy.is_finite():
            raise ValueError('Invalid measurement or energy')
        if (state, seed) in by_ligand[ligand] or measurements.setdefault(ligand, nm) != nm:
            raise ValueError('Duplicate observation or inconsistent label for ' + ligand)
        by_ligand[ligand][state, seed] = energy
    blocks = set(next(iter(by_ligand.values())))
    states = sorted({state for state, _ in blocks})
    seeds = sorted({seed for _, seed in blocks})
    if len(states) != 2 or blocks != set(itertools.product(states, seeds)):
        raise ValueError('Expected a complete two-state design')
    if any(set(data) != blocks for data in by_ligand.values()):
        raise ValueError('Unequal receptor or replicate coverage')
    energy = {l: tuple(float(statistics.median(data[s, seed] for seed in seeds)) for s in states)
              for l, data in by_ligand.items()}
    origin = min(e for pair in energy.values() for e in pair)
    ligands = []
    for l, values in sorted(energy.items()):
        k = tuple(Fraction(math.exp(-(e - origin) / RT)) for e in values)
        if min(k) == 0:
            raise ValueError('Association proxy underflows')
        ligands.append(Ligand(l, measurements[l], values, k))
    return states, seeds, ligands


def sign_credit(value):
    return (1 + (value > 0) - (value < 0)) / 2


def crossings(edges):
    """Exact sweep in rational arithmetic after converting energies to binary64 activities.

    Each edge has an affine activity difference with values d0, d1 at p=0,1.
    Include isolated ties at the endpoints. Concurrent interior crossings split
    their scores halfway across the jump, even for opposing edge directions.
    """
    jumps = defaultdict(int)
    left = right = inside = separately = 0.
    for d0, d1 in edges:
        c0, c1 = sign_credit(d0), sign_credit(d1)
        left += c0
        right += c1
        separately += max(c0, c1)
        inside += c0 if d0 else c1
        if d0 * d1 < 0:
            jumps[d0 / (d0 - d1)] += 1 if d1 > 0 else -1
    best = max(left, right, inside)
    for point in sorted(jumps):
        inside += jumps[point]
        best = max(best, inside)
    return separately, best


def analyze(key, records):
    states, seeds, ligands = prepare(records)
    info = dict(zip(('target','assay','endpoint','reference'), key))
    edges, details, floors, errors = [], [], [], []
    observed = [float(l.nm) for l in ligands]
    mixture = [-RT * math.log(float(sum(l.activity) / 2)) for l in ligands]
    lower_constraint, upper_constraint = None, None
    for a, b in itertools.combinations(sorted(ligands, key=lambda l: (l.nm, l.name)), 2):
        if b.nm < 3*a.nm:
            continue
        exp_gap = RT * math.log(float(b.nm / a.nm))
        gaps = tuple(y-x for x,y in zip(a.energy, b.energy))
        low, high = min(gaps), max(gaps)
        gap = RT * math.log(float(sum(a.activity) / sum(b.activity)))
        # The original per-pair implementation uses independent log-domain energy aggregation.
        existing = pair_sensitivity(a.energy, b.energy)['gap_equal_population']
        if abs(gap - existing) > 1e-10:
            raise ArithmeticError('Failed per-pair regression')
        d1, d0 = (x-y for x,y in zip(a.activity,b.activity))
        edges.append((d0, d1))
        floor = max(low-exp_gap, exp_gap-high, 0.)
        floors.append(floor)
        errors.append(gap-exp_gap)
        bound = ''
        direction = 'constant'
        if d0*d1 < 0:
            point = d0/(d0-d1)
            bound = float(point)
            if d1 > 0:
                direction = 'above'
                lower_constraint = max(lower_constraint or (point,a.name,b.name), (point,a.name,b.name))
            else:
                direction = 'below'
                upper_constraint = min(upper_constraint or (point,a.name,b.name), (point,a.name,b.name))
        details.append(dict(**info,strong=a.name,weak=b.name,fold=float(b.nm/a.nm),experimental_gap=exp_gap,
                            gap_at_p0=gaps[1],gap_at_p1=gaps[0],gap_equal=gap,
                            equal_credit=sign_credit(d0+d1),best_separate_credit=max(sign_credit(d0),sign_credit(d1)),
                            interior_correct=direction,crossover=bound,min_absolute_error=floor))
    if not edges:
        raise ValueError(f'No eligible pair: {key}')
    separate, common = crossings(edges)
    n = len(edges)
    conflict = (lower_constraint and upper_constraint and lower_constraint[0] >= upper_constraint[0])
    row = dict(**info,state_at_p1=states[0],state_at_p0=states[1],seeds=len(seeds),ligands=len(ligands),pairs=n,
               equal_credit=sum(d['equal_credit'] for d in details),separate_ceiling=separate,shared_ceiling=common,
               fixed_order_loss=n-separate,common_prior_loss=separate-common,
               rho_equal=spearmanr(observed,mixture).statistic,tau_b_equal=kendalltau(observed,mixture).statistic,
               gap_mae_equal=statistics.mean(map(abs,errors)),gap_mae_floor=statistics.mean(floors),
               gap_rmse_equal=math.sqrt(statistics.mean(x*x for x in errors)),
               gap_rmse_floor=math.sqrt(statistics.mean(x*x for x in floors)),gap_bias_equal=statistics.mean(errors),
               witness_above=('>'.join(lower_constraint[1:]) if conflict else ''),
               witness_below=('>'.join(upper_constraint[1:]) if conflict else ''))
    return row, details


def main(out):
    panels, paths = collect()
    summary, pairs = [], []
    for key, records in sorted(panels.items()):
        row, edges = analyze(key, records)
        summary.append(row)
        pairs.extend(edges)
    out.mkdir(parents=True, exist_ok=True)
    write(out/'summary.csv',summary)
    write(out/'pairs.csv',pairs)
    hashes = {str(p.relative_to(BENCH.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (out/'inputs.json').write_text(json.dumps(hashes,indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'results')
    main(parser.parse_args().output)
