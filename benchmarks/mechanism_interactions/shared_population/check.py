"""Independent two-state mixture replay, edge bounds, witness and precision checks."""
import argparse
import csv
import hashlib
import itertools
import json
import lzma
import math
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from statistics import mean, median

from scipy.stats import kendalltau, spearmanr

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
RT = Decimal('0.00198720425864083') * Decimal('298.15')
KEY = ('target', 'assay', 'endpoint', 'reference')


def read(path):
    opener = lzma.open if path.suffix == '.xz' else open
    with opener(path, 'rt') as f:
        return list(csv.DictReader(f))


def observations():
    directly_printed = {}
    raw = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    labels = {}
    kd = BENCH / 'kinase_state_affinity'
    panel = {(r['target'], r['cid']): r['kd_nm'] for r in read(kd / 'panel.csv')}

    def add(key, ligand, receptor, seed, nm, score, native_precision):
        replicates = raw[key][ligand][receptor]
        if seed in replicates or labels.setdefault((key, ligand), Decimal(nm)) != Decimal(nm):
            raise ValueError('Duplicate seed or inconsistent endpoint within ligand and assay')
        replicates[seed] = Decimal(score)
        directly_printed[key] = native_precision

    for table in ('results/cells.csv', 'independent-seeds/results/cells.csv'):
        for r in read(kd / table):
            target = r['target']
            for reference, field in [('prepared', 'energy'), ('pool', 'native')]:
                key = (target, f'Davis-{target}', 'Kd', reference)
                add(key, r['ligand'], r['state'], r['replicate'], panel[target, r['ligand']],
                    r[field], field == 'native')

    path = BENCH / 'frozen_core/strain_rotor_interaction/first_poses.csv.xz'
    included = {'fxa_repaired', 'native_12k', 'transfer_archived', 'transfer_active_torsions', 'parp_spr'}
    for r in read(path):
        if r['panel'] in included:
            reference = 'active_rotor_control' if r['panel'] == 'transfer_active_torsions' else 'pool'
            key = (r['target'], r['assay'], r['endpoint'], reference)
            add(key, r['compound'], r['receptor'], r['seed'], r['measured_nM'],
                r['native'], r['panel'] != 'transfer_active_torsions')
    return raw, labels, directly_printed


def checked_medians(panel, states):
    medians, seed_set = {}, None
    for ligand, receptor_data in panel.items():
        assert set(states) == set(receptor_data)
        for scores in receptor_data.values():
            if seed_set is None:
                seed_set = set(scores)
            assert seed_set == set(scores)
        medians[ligand] = tuple(median(receptor_data[s].values()) for s in states)
    return medians, len(seed_set)


def credit(value):
    return (int(value > 0) + int(value >= 0)) / 2


def activities(energies):
    # Direct evaluation in high precision from exact decimal input medians.
    return {l: tuple(Fraction((-e / RT).exp()) for e in pair) for l, pair in energies.items()}


def switches(K, edges, strong_factor=Fraction(1)):
    points = {Fraction(0), Fraction(1)}
    for a, b in edges:
        d0, d1 = [strong_factor*x-y for x, y in zip(K[a], K[b])]
        if d0 != d1:
            p = d0/(d0-d1)
            if 0 < p < 1:
                points.add(p)
    cuts = sorted(points)
    return cuts + [(lo+hi)/2 for lo, hi in zip(cuts, cuts[1:])]


def ceiling(K, edges, strong_factor=Fraction(1)):
    # Directly construct each ligand's mixture, evaluate credit at each cut and
    # within every segment. This does not count directional jumps.
    best = 0
    for p in switches(K, edges, strong_factor):
        mixture = {l: (1-p)*k0+p*k1 for l, (k0,k1) in K.items()}
        if strong_factor == 1:
            rank = {value: i for i, value in enumerate(sorted(set(mixture.values())))}
            ranks = {l: rank[v] for l, v in mixture.items()}
            n = sum(credit(ranks[a] - ranks[b]) for a,b in edges)
        else:
            strong = {l: value*strong_factor for l, value in mixture.items()}
            n = sum(0.5 * (int(strong[a] > mixture[b]) + int(strong[a] >= mixture[b])) for a,b in edges)
        best = max(best,n)
    return best


def close(actual, expected):
    if not math.isclose(float(actual), float(expected), rel_tol=1e-10, abs_tol=1e-10):
        raise ValueError(f'Disagreement: {actual!r}, {expected!r}')


def validate_panel(row, pair_rows, key, raw, labels, printed):
    states = sorted({s for receptor_map in raw.values() for s in receptor_map}, reverse=True)
    assert states == [row['state_at_p0'], row['state_at_p1']]
    energy, nseed = checked_medians(raw, states)
    assert len(energy) == int(row['ligands']) and nseed == int(row['seeds'])
    nM = {l: labels[(key,l)] for l in energy}
    ordered = sorted(energy, key=nM.get)
    edges = [(a,b) for a,b in itertools.combinations(ordered,2) if nM[b] >= 3*nM[a]]
    assert len(edges) == len(pair_rows) == int(row['pairs'])
    by_edge = {(r['strong'],r['weak']):r for r in pair_rows}
    assert set(edges) == by_edge.keys()
    K = activities(energy)
    equal_mix = {l:(x+y)/2 for l,(x,y) in K.items()}
    observed = [nM[l] for l in ordered]
    strength = [-float(equal_mix[l]) for l in ordered]  # monotone with energy
    close(spearmanr(observed,strength).statistic, row['rho_equal'])
    close(kendalltau(observed,strength).statistic, row['tau_b_equal'])
    pure, equal, errors, floors = [], [], [], []
    lower, upper = [], []
    for a,b in edges:
        r = by_edge[a,b]
        expected = RT * (nM[b]/nM[a]).ln()
        gaps = tuple(float(v-u) for u,v in zip(energy[a],energy[b]))
        mixed = float(RT) * math.log(float(equal_mix[a]/equal_mix[b]))
        floor = max(min(gaps) - float(expected), float(expected) - max(gaps), 0.)
        close(float(nM[b]/nM[a]),r['fold'])
        close(expected,r['experimental_gap'])
        close(gaps[0],r['gap_at_p0'])
        close(gaps[1],r['gap_at_p1'])
        close(mixed,r['gap_equal'])
        close(floor,r['min_absolute_error'])
        pure.append(max(credit(gaps[0]),credit(gaps[1])))
        equal.append(credit(equal_mix[a]-equal_mix[b]))
        close(pure[-1],r['best_separate_credit'])
        close(equal[-1],r['equal_credit'])
        floors.append(floor)
        errors.append(mixed-float(expected))
        if min(gaps) < 0 < max(gaps):
            point = next(p for p in switches(K,[(a,b)]) if 0 < p < 1)
            (lower if gaps[1] > 0 else upper).append((point, f'{a}>{b}'))
            assert r['interior_correct'] == ('above' if gaps[1] > 0 else 'below')
            close(point,r['crossover'])
        else:
            assert r['crossover'] == '' and r['interior_correct'] == 'constant'
    best = ceiling(K,edges)
    separate = sum(pure)
    for field,value in [('equal_credit',sum(equal)),('separate_ceiling',separate),('shared_ceiling',best),
                        ('fixed_order_loss',len(edges)-separate),('common_prior_loss',separate-best),
                        ('gap_mae_equal',mean(map(abs,errors))),('gap_mae_floor',mean(floors)),
                        ('gap_rmse_equal',math.sqrt(mean(x*x for x in errors))),
                        ('gap_rmse_floor',math.sqrt(mean(x*x for x in floors))),('gap_bias_equal',mean(errors))]:
        close(value,row[field])
    if lower and upper and max(lower)[0] >= min(upper)[0]:
        L = dict((pair,p) for p,pair in lower)[row['witness_above']]
        U = dict((pair,p) for p,pair in upper)[row['witness_below']]
        assert L == max(lower)[0] and U == min(upper)[0] and L >= U
    else:
        assert not row['witness_above'] and not row['witness_below']
    # Each pair is optimistically allowed its OWN +/- energy shifts. Sharing a
    # ligand across different edges need not permit all these shifts together.
    factor = Fraction((Decimal('0.001')/RT).exp())
    relaxed = ceiling(K,edges,factor)
    separate_relaxed = sum(max(credit(factor*x-y) for x,y in zip(K[a],K[b])) for a,b in edges)
    assert best <= relaxed <= separate_relaxed
    return dict(**dict(zip(KEY,key)), pairs=len(edges),
                uncertainty='print_half_step' if printed else 'hypothetical_energy_shift',
                half_width_per_ligand='0.0005',shared_credit_nominal=best,
                shared_credit_optimistic=relaxed,separate_credit_optimistic=separate_relaxed,
                unresolved_credit_at_least=len(edges)-relaxed)


def check_controls():
    # Only the synthetic controls call the implementation under review. All
    # experimental panel statistics above are reconstructed independently.
    from replay import crossings, prepare

    # The chosen population 1/4 gives integer association constants for this
    # constructed assay. Wider-than-threefold pairs have the known direction.
    K = {l:tuple(map(Fraction,k)) for l,k in {'A':(20,8),'B':(4,1),'C':(1,2)}.items()}
    truth = {l:Fraction(3,4)*x+Fraction(1,4)*y for l,(x,y) in K.items()}
    edges = [(a,b) for a,b in itertools.permutations(K,2) if truth[a]>=3*truth[b]]
    assert len(edges) >= 2 and ceiling(K,edges) == len(edges)
    # Competing directions with a simultaneous crossing, identical lines,
    # endpoint-only ties and a dependence impossible at any finite population.
    ties = {l:tuple(map(Fraction,k)) for l,k in
            {'A':(1,3),'B':(2,2),'C':(2,2),'D':(1,1),'E':(3,3)}.items()}
    tests = [([('A','B'),('C','A')],1),([('B','C')],.5),([('A','E')],.5),
             ([('D','A')],.5),([('D','E')],0)]
    for edges,expected in tests:
        assert ceiling(ties,edges) == expected
        swapped = {l:p[::-1] for l,p in ties.items()}
        scaled = {l:(a*7,b*7) for l,(a,b) in ties.items()}  # overall gauge
        # A reference shift common to all ligands within each state induces a
        # one-to-one reparametrization of p. The feasible union is invariant.
        bystate = {l:(a*11,b*3) for l,(a,b) in ties.items()}
        for states in (ties, swapped, scaled, bystate):
            assert ceiling(states,edges) == expected
            deltas = [tuple(x-y for x,y in zip(states[a],states[b])) for a,b in edges]
            assert crossings(deltas)[1] == expected

    good = [(l,s,seed,label,score) for l,label,score in [('A','1','-2'),('B','9','-1')]
            for s in ('R1','R2') for seed in ('1','2')]
    assert len(prepare(good)[2]) == 2
    invalid = [good + [good[0]], good[:-1],
               [(*good[0][:3],'77',good[0][-1])] + good[1:]]
    for bad in invalid:
        try:
            prepare(bad)
        except ValueError:
            continue
        raise AssertionError('Duplicate, incomplete block or contradicted label was accepted')


def main(results, precision_out):
    check_controls()
    raw, labels, printed = observations()
    all_pairs = defaultdict(list)
    for r in read(results / 'pairs.csv'):
        all_pairs[tuple(r[k] for k in KEY)].append(r)
    summaries = read(results / 'summary.csv')
    by_key = {tuple(r[k] for k in KEY):r for r in summaries}
    eligible = {k for k,panel in raw.items() if len({s for l in panel.values() for s in l}) == 2}
    assert eligible == by_key.keys() == all_pairs.keys() and len(by_key) == len(summaries)
    for p,digest in json.loads((results / 'inputs.json').read_text()).items():
        assert hashlib.sha256((BENCH.parent/p).read_bytes()).hexdigest() == digest
    precision = [validate_panel(by_key[k],all_pairs[k],k,raw[k],labels,printed[k]) for k in sorted(eligible)]
    if precision_out:
        with precision_out.open('w', newline='') as f:
            w = csv.DictWriter(f,precision[0],lineterminator='\n')
            w.writeheader()
            w.writerows(precision)
    print(f'Independent direct mixture: {len(precision)} panels, {sum(len(r) for r in all_pairs.values())} pair rows.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,default=HERE/'results')
    parser.add_argument('--precision-output',type=Path)
    args = parser.parse_args()
    with localcontext() as ctx:
        ctx.prec = 50
        main(args.results, args.precision_output)
