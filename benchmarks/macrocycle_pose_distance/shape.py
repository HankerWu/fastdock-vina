"""Diagnose residual rigid placement and intraligand shape on every released pose."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "macrocycle_proposal" / "external"))
from external_score import poses, write


def fitted_rmsd(centered, target):
    u, _, vt = np.linalg.svd(centered.T @ target)
    flip = np.sign(np.linalg.det(u @ vt))
    r = (u * [1, 1, flip]) @ vt
    return np.sqrt(np.sum((centered @ r - target)**2) / len(centered))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs",type=Path)
    p.add_argument("results",type=Path)
    a=p.parse_args()
    rows=[]
    for path in sorted((a.results/"raw").glob("*.pdbqt")):
        pdb = path.name[:4]
        meta=json.loads((a.inputs/pdb/"preparation.json").read_text())
        refs=np.asarray(meta["mapped"])
        ref0=refs - refs.mean(axis=1,keepdims=True)
        for rank,pose in enumerate(poses(path),1):
            xyz=np.asarray([pose[i][1] for i in meta["serial"]])
            centroid=float(np.linalg.norm(xyz.mean(axis=0)-refs[0].mean(axis=0)))
            x0=xyz-xyz.mean(axis=0)
            receptor=np.sqrt(np.sum((xyz-refs)**2,axis=(1,2))/len(xyz)).min()
            recentered=np.sqrt(np.sum((x0-ref0)**2,axis=(1,2))/len(xyz)).min()
            shape=min(fitted_rmsd(x0, t) for t in ref0)
            assert np.isclose(receptor**2,centroid**2+recentered**2) and shape <= recentered + 1e-9
            rows.append(dict(key=path.stem,rank=rank,rmsd_A=receptor,centroid_error_A=centroid,
                             recentered_rmsd_A=recentered,rigid_fitted_rmsd_A=shape))
    write(a.results/"shape.csv",rows)
