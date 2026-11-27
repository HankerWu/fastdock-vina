"""Post hoc lower bounds on the best achievable distance, independent of the torsion solver."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from manifold import tree
from verify import BENCH, Ligand, REC, active_pools, members, modes


def block_residual(x, y, weights):
    """Relax a rigid block to any proper rotation and translation."""
    source = x - np.average(x, weights=weights, axis=0)
    dest = y - np.average(y, weights=weights, axis=0)
    u, _, vt = np.linalg.svd((weights[:, None]*source).T @ dest)
    u[:, -1] *= np.linalg.det(u @ vt)
    residual = np.square(source @ (u @ vt)-dest).sum(axis=1)
    return residual @ weights


def main(out):
    path = BENCH / "cubic_grid/human_panel/affinity-pools.tar.xz"
    prep = members(path, ["affinity/5.pdbqt","affinity/6.pdbqt"])
    five,six = [Ligand(prep[f"affinity/{i}.pdbqt"]) for i in (5,6)]
    _, tree_groups, blocks, _ = tree(prep["affinity/5.pdbqt"],five)
    shares = 1/np.bincount(np.concatenate(blocks))
    mappings = six.graph.GetSubstructMatches(five.graph, useChirality=True, uniquify=False)
    frame = pd.read_csv(out/"manifold_poses.csv")
    data = np.load(out/"manifold_xyz.npz")
    assert data["pdbqt_serials"].tolist() == five.order
    assert len(frame)==len(data["xyz"])==234
    rows=[]
    for panel, arm, seed, pools in active_pools():
        for rec in REC:
            for rank, lines in enumerate(modes(pools[rec,6]),1):
                row = frame.iloc[len(rows)]
                assert (row.panel,row.arm,row.seed,row.receiver,row["rank"]) == (panel,arm,seed,rec,rank)
                xyz = data["xyz"][len(rows)]
                target = six.coords(lines)
                ub = min(np.linalg.norm(xyz-target[list(m)],"fro")/np.sqrt(len(m)) for m in mappings)
                assert abs(ub-row.feasible_rmsd) < 1e-9
                # Each atom's repeated weight across blocks is exactly 1.
                # Separate rigid fits can only improve on any one feasible tree.
                costs = np.array([[block_residual(five.reference[group],target[list(m)][group],shares[group])
                                   for group in blocks] for m in mappings])
                best = int(np.argmin(costs.sum(axis=1)))
                total = costs[best].sum()
                lb = np.sqrt(total/len(five.order))
                # The six heavy atoms of the pipecolate ring form 5's ROOT.
                root = tree_groups[0]
                ring_rmsd = np.sqrt(block_residual(five.reference[root], target[list(mappings[best])][root],
                                                   np.ones(len(root)))/len(root))
                if lb>ub+1e-9:
                    raise ArithmeticError("Bound exceeds a measured feasible witness")
                rows.append(dict(panel=panel,arm=arm,seed=seed,receiver=rec,rank=rank,
                                 lower_rmsd=lb, upper_rmsd=ub,gap=ub-lb,
                                 piperidine_ring_only_rmsd=ring_rmsd,
                                 root_and_flanking_axis_fraction=costs[best,0]/total))
    result=pd.DataFrame(rows)
    result.to_csv(out/"manifold_bounds.csv",index=False)
    details=dict(poses=len(rows), groups=len(blocks), max_gap=float(result.gap.max()),
                 lower_range=[float(result.lower_rmsd.min()),float(result.lower_rmsd.max())],
                 median_lower=float(result.lower_rmsd.median()), median_gap=float(result.gap.median()))
    print(json.dumps(details,indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out",type=Path)
    main(parser.parse_args().out)
