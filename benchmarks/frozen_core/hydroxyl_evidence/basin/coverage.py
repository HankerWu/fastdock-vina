"""Rigid best-fit lower bound on crystal RMSD of every emitted blind internal conformation."""
import argparse
import json
from pathlib import Path

import numpy as np

from swap import coords, rmsd, transform, write


def main(work):
    source=work/"source"
    rows=[]
    for pdb in ("2YIY","4I1R"):
        meta=json.loads((source/pdb/"preparation.json").read_text())
        refs=np.array(meta["mapped"])
        for pool in sorted((source/"runs").glob(pdb+"_*_pro.pdbqt")):
            _,cap,seed,_=pool.stem.split("_")
            for rank,block in enumerate((b for b in pool.read_text().split("ENDMDL") if b.strip()),1):
                points=coords(block)
                xyz=np.array([points[i] for i in meta["serial"]])
                errors=[]
                for ref in refs:
                    r,t=transform(xyz,ref)
                    errors.append(float(np.sqrt(((xyz@r+t-ref)**2).sum(axis=1).mean())))
                rows.append(dict(pdb=pdb,cap=cap,seed=seed,rank=rank,blind_rmsd=rmsd(xyz,refs),
                                 minimum_over_rigid_motion=min(errors),centroid_distance=np.linalg.norm(xyz.mean(0)-refs[0].mean(0))))
    write(work/"internal_coverage.csv",rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    main(p.parse_args().work)
