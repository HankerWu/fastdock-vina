"""Evaluate the complete legal and in-box domains for compounds with deposited coordinates."""
import argparse
from pathlib import Path
import sys
import numpy as np

from trial import HERE, LIVE, native

sys.path.insert(0,str(HERE.parent/"engine"))
from coverage_gate import candidates
from pocket_triangle import Pocket
sys.path.insert(0,str(native.ROOT/"benchmarks/cavity_field_rank/fkbp"))
from frozen_assay import old
from check_assay import rms


def evaluate(work,out):
    # check_panel extracts these receiver/CCD records from the hash-verified deposition archive.
    refs={r:old.reference(work,r,c) for r,c in zip(native.RECEPTORS,("SB3","SBX"),strict=True)}
    rows=[]
    for rec in native.RECEPTORS:
        with (work/"maps"/rec/"receptor.C.map").open() as f: h=[next(f).split() for _ in range(6)]
        pocket=Pocket.from_pdbqt(work/rec/"receptor.pdbqt", np.array(h[5][1:],float),
                                 float(h[3][1])*np.array(h[4][1:],int))
        for n,crystal_rec in ((8,"1FKG"),(9,"1FKH")):
            mol,chem_order,_=old.molecule((work/"affinity"/f"{n}.pdbqt").read_text())
            ref=old.crystal_positions(refs[crystal_rec],refs[rec],mol)
            stem=f"{rec}_{n}"
            source=(work/"pools/edge"/stem).with_suffix(".source").read_text().splitlines()
            serial=np.array(source[0][1:].split(),int)
            order=np.argsort(serial)
            x=np.array(source[1].split()[3:],float).reshape(-1,3)[order]
            serial=serial[order]
            chem=np.array([serial.tolist().index(i) for i in chem_order])
            pool=candidates(pocket,x,serial)
            errors=[]
            for ids in np.array_split(pool["triangle_id"],150):
                xyz=(x-x.mean(0)) @ pool["matched_rotation"][ids].transpose(0,2,1) + pool["center"][ids,None]
                errors.extend(rms(xyz[:,chem],ref))
            err=np.array(errors)
            eligible=pool["eligible"]
            for arm in LIVE:
                with np.load((work/"pools"/arm/stem).with_suffix(".npz")) as b:
                    ids=b["triangle_id"]
                    assert np.array_equal(pool["site_triangles"][ids],b["site_triangles"])
                    assert np.allclose(pool["center"][ids],b["center"],atol=1e-12,rtol=0)
                    assert np.allclose(pool["matched_rotation"][ids],b["rotation"],atol=1e-12,rtol=0)
                    ranked=ids[np.argsort(b["true_score"],kind="stable")]
                    if arm!="edge": assert eligible[ids].all()
                    near=np.flatnonzero(err[ranked]<=2)
                    rows.append(dict(receiver=rec,compound=n,arm=arm,all=len(err),feasible=eligible.sum(),
                        all_min_A=err.min(),feasible_min_A=err[eligible].min(),bank_min_A=err[ids].min(),
                        all_le2=(err<=2).sum(),feasible_le2=(err[eligible]<=2).sum(),bank_le2=len(near),
                        best_field_near_rank=near[0]+1 if len(near) else "",winner_A=err[ranked[0]]))
    native.table(out,rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    evaluate(a.work.resolve(),a.out.resolve())
