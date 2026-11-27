"""Verify both Vina-reported RMSD bounds directly from the released Cartesian coordinates."""
import argparse
from pathlib import Path
import sys

import numpy as np
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "macrocycle_proposal" / "external"))
from external_score import poses, write


def element(t):
    if t.startswith("CG") or t in ("A", "C"):
        return "C"
    if t.startswith("G") or t == "W":
        return "dummy"
    if t in ("NA", "OA", "SA"):
        return t[0]
    return t


def bounds(pose, ref, include_markers):
    ids = [i for i, (t, _) in ref.items() if t not in ("HD", "H") and (include_markers or not t.startswith("G"))]
    xyz = np.array([pose[i][1] for i in ids])
    xyz0 = np.array([ref[i][1] for i in ids])
    d2 = cdist(xyz,xyz0,metric="sqeuclidean")
    el = np.array([element(ref[i][0]) for i in ids])
    d2[el[:,None] != el[None,:]] = np.inf
    lb = np.sqrt(max(d2.min(axis=1).mean(), d2.min(axis=0).mean()))
    return lb, np.sqrt(np.diag(d2).mean())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results", type=Path)
    a = p.parse_args()
    rows = []
    for path in sorted((a.results/"raw").glob("*.pdbqt")):
        block = list(poses(path))
        scores = [l.split()[-2:] for l in path.read_text().splitlines() if l.startswith("REMARK VINA RESULT:")]
        for rank,(pose,score) in enumerate(zip(block, scores),1):
            old = bounds(pose,block[0],True)
            fixed = bounds(pose,block[0],False)
            official = np.array(score,dtype=float)
            expected = fixed if path.stem.endswith("_fixed") else old
            # Coordinates and output distances both rounded to 0.001 Angstrom.
            assert np.max(np.abs(expected-official)) < 0.00224
            rows.append(dict(key=path.stem,rank=rank,printed_lb=score[0],printed_ub=score[1],
                             old_lb=old[0],old_ub=old[1],fixed_lb=fixed[0],fixed_ub=fixed[1]))
    write(a.results/"reported.csv",rows)
