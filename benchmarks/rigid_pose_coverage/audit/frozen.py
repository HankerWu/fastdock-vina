"""Post-prediction audit of all frozen historical whole-ligand correspondences."""
import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
sys.path[:0] = [str(HERE.parent / 'engine'), str(BENCH / 'fragment_anchor/audit')]
from coverage_gate import candidates, embedding, select
from edge_audit import Pocket, equivalent_targets, frames, raw_source, table


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def distances_to_crystal(x, bank, targets):
    centered = x - x.mean(0)
    result = np.empty(len(bank['triangle_id']))
    for start in range(0, len(result), 128):
        region = slice(start, start + 128)
        coords = centered @ bank['matched_rotation'][region].transpose(0, 2, 1)
        coords += bank['center'][region, None, :]
        result[region] = np.sqrt(np.minimum.reduce(
            [np.square(coords - target).sum(axis=2).mean(axis=1) for target in targets]))
    return result


def audit(p, work, out):
    t0, c0 = time.perf_counter(), time.process_time()
    out.mkdir(parents=True, exist_ok=True)
    data = frames(work)
    manifest = json.loads((p/'inputs.json').read_text())['cells']
    historical = {r['case']: r for r in csv.DictReader(
        (BENCH/'fragment_anchor/audit/results/gates_pose.csv').open()) if r['arm']=='whole'}
    rows, witnesses, banks = [], [], {}
    for case, old in sorted(historical.items()):
        name, seed = case.rsplit('_', 1)
        spec = manifest[name]
        assert all(digest(p/path) == expected for path, expected in spec['sha256'].items())
        source = p/'dump'/(case+'_whole')/'source.xyz'
        ids, x = raw_source(source)
        pocket = Pocket.from_pdbqt(p/spec['receptor'], spec['center'], spec['size'])
        bank = candidates(pocket, x, ids)
        allowed = np.flatnonzero(bank['eligible'])
        feats = embedding(x, bank['center'], bank['matched_rotation'])
        # Only this analysis reads the label. All prediction and selection inputs are now fixed.
        order, refs = equivalent_targets(name, p/spec['ligand'], data)
        distances = distances_to_crystal(x, bank, refs[:, [order.index(i) for i in ids]])
        np.testing.assert_allclose(distances.min(), float(old['before_min_A']), atol=1e-10, rtol=0)
        near = np.flatnonzero(distances <= 2.)
        assert len(near) == int(old['before_le2'])
        for arm in ('edge', 'feasible_edge', 'coverage', 'uniform'):
            selected, radius2 = select(bank, x, int(seed), arm)
            if arm == 'edge':
                np.testing.assert_allclose(distances[selected].min(),float(old['after_min_A']),atol=1e-10,rtol=0)
            nearest = cKDTree(feats[selected])
            d, _ = nearest.query(feats[allowed], workers=1)
            if np.isfinite(radius2):
                np.testing.assert_allclose(np.sqrt(radius2), d.max(), atol=1e-10, rtol=0)
            oracle = distances[selected].min()
            assert oracle <= distances[allowed].min() + d.max() + 1e-10
            rows.append(dict(case=case,arm=arm,source_sha256=digest(source),all_count=len(distances),
                feasible_count=len(allowed),selected_count=len(selected),selected_outside=int(sum(~bank['eligible'][selected])),
                ideal_rigid_A=old['ideal_rigid_A'],all_min_A=distances.min(),all_near=len(near),
                feasible_min_A=distances[allowed].min(),feasible_near=int(sum(distances[allowed] <= 2)),
                selected_min_A=oracle,selected_near=int(sum(distances[selected] <= 2)),
                nearest_median_A=np.median(d),nearest_P90_A=np.quantile(d,.9),
                nearest_P95_A=np.quantile(d,.95),radius_A=d.max(),
                selected_best_triangle_id=int(selected[np.argmin(distances[selected])])))
            banks[case+'_'+arm] = selected.astype(np.int32)
            for i in near:
                delta, nn = nearest.query(feats[i])
                j = selected[nn]
                witnesses.append(dict(case=case,arm=arm,near_id=i,feasible=int(bank['eligible'][i]),
                    crystal_A=distances[i],selected=int(i in selected),nearest_selected_id=j,
                    labeled_distance_A=delta,neighbor_crystal_A=distances[j],edge_residual_A2=bank['edge_residual_A2'][i]))
    table(out/'historical.csv', rows)
    table(out/'near_neighbors.csv', witnesses)
    np.savez_compressed(out/'selected_ids.npz', **banks)
    (out/'cost.json').write_text(json.dumps(dict(wall_s=time.perf_counter()-t0,cpu_s=time.process_time()-c0,
        n_source=len(historical),n_arms=4,includes='offline enumeration, selection, crystal and distance audits'),indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('parent_pose', 'frame_work', 'results'):
        parser.add_argument(arg,type=Path)
    args = parser.parse_args()
    audit(args.parent_pose,args.frame_work,args.results)
