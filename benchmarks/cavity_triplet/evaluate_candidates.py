"""Reconstruct every blind pose after prediction, using two independent experimental references."""

import argparse
import csv
import gzip
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
from rdkit.Chem import rdMolAlign
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

from generate import Pocket, candidate_coords, write

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "receptor_linkage/aurora/internal_loops/rigid_coupled"))
import study


def check(source_table, candidates, aurora, brd, dest):
    dest.mkdir(exist_ok=True,parents=True)
    with tempfile.TemporaryDirectory() as tmp:
        primary = study.scan.config(aurora,brd,Path(tmp))
        independent = study.audit.inputs(aurora,brd,Path(tmp))
    cells, maximum, pockets = [], dict(crystal_A=0., frame_A=0., frame_score=0., direct_score=0.), {}
    # Rigid change of the source coordinates only; fixed before evaluation.
    move = Rotation.from_rotvec([0.21,0.77,-1.53])
    offset = np.array([15.25,-3.,2.])
    detail_path = dest / "candidate_metrics.csv.gz"
    with gzip.open(detail_path, "wt",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["source","arm","index","score","rmsd_A","coarse_rank"])
        for row in csv.DictReader(source_table.open()):
            source, drug, receptor = (row[k] for k in ("source","drug","receptor"))
            with np.load(row["npz"]) as parent:
                x, serial = parent["xyz"], parent["serials"].tolist()
            center = tuple(float(row["center_"+a]) for a in "xyz")
            k = receptor, center, float(row["size"])
            if k not in pockets:
                pockets[k] = Pocket.from_pdbqt(*k)
            pocket = pockets[k]
            refs = primary[drug]["refs"][:, [primary[drug]["serials"].index(s) for s in serial]]
            mol, order, ref = independent[drug][-3:]
            independent_ref = ref.GetConformer().GetPositions()
            maps = ref.GetSubstructMatches(mol, uniquify=False, maxMatches=100000)
            if not maps or len(maps)>=100000:
                raise ValueError("reference mapping missing or saturated")
            independent_refs = independent_ref[np.asarray(maps)]
            reorder = [serial.index(s) for s in order]
            other_x = move.apply(x) + offset
            transformed = pocket.place(other_x,int(row["seed"]))
            with np.load(candidates / (source+".npz")) as res:
                for arm in ("matched","sham"):
                    xyz = candidate_coords(x,res,arm)
                    rmsd = np.sqrt(np.square(xyz[:,None,:,:]-refs).sum(axis=-1).mean(axis=-1)).min(axis=1)
                    direct_rmsd = np.sqrt(np.square(xyz[:,None,reorder,:]-independent_refs).sum(axis=-1).mean(axis=-1)).min(axis=1)
                    maximum["crystal_A"] = max(maximum["crystal_A"],abs(rmsd-direct_rmsd).max())
                    maximum["frame_A"] = max(maximum["frame_A"],abs(xyz-candidate_coords(other_x,transformed,arm)).max())
                    score = res[arm+"_score"]
                    maximum["frame_score"] = max(maximum["frame_score"],abs(score-transformed[arm+"_score"]).max())
                    priority = np.argsort(np.round(score,10),kind="stable")
                    rank = np.argsort(priority) + 1
                    chosen = int(res[arm+"_best"])
                    assert chosen==transformed[arm+"_best"]==priority[0]
                    assert np.array_equal(res["site_triangles"],transformed["site_triangles"])
                    # Verify the score of the chosen and label-closest candidates by all-pairs distance.
                    for j in {chosen,int(rmsd.argmin())}:
                        u = cdist(xyz[j],pocket.tree.data).min(axis=1) - 3.8
                        outside = np.maximum(np.abs(xyz[j]-center)-float(row["size"])/2,0.)
                        independent_score = sum(max(0.,-s)**2-np.exp(-s**2) for s in u) + (outside**2).sum()
                        maximum["direct_score"] = max(maximum["direct_score"],abs(score[j]-independent_score))
                        for i, pos in enumerate(xyz[j,reorder]):
                            mol.GetConformer().SetAtomPosition(i,pos)
                        direct = rdMolAlign.CalcRMS(mol,ref,symmetrizeConjugatedTerminalGroups=False)
                        maximum["crystal_A"] = max(maximum["crystal_A"],abs(direct-rmsd[j]))
                    writer.writerows((source,arm,i,float(score[i]),float(rmsd[i]),int(rank[i])) for i in range(len(rmsd)))
                    near = rmsd<=2
                    cells.append(dict(source=source, drug=drug, seed=row["seed"],trajectory=row["trajectory"],step=row["step"],arm=arm,
                                      candidate_count=len(rmsd),generator_best_A=rmsd.min(),
                                      candidates_le2=int(near.sum()), candidates_le3=int((rmsd<=3).sum()),
                                      best_candidate_coarse_rank=rank[rmsd.argmin()],first_le2_coarse_rank=rank[near].min() if near.any() else None,
                                      selected_A=rmsd[chosen],selected_index=chosen,selected_score=score[chosen]))
    if max(maximum.values())>=1e-7:
        raise ValueError(str(maximum))
    write(dest / "candidate_stages.csv",cells)
    (dest / "candidate_audit.json").write_text(json.dumps(dict(max_errors=maximum, source_events=len(cells)//2,
        candidates=sum(r["candidate_count"] for r in cells), held_parent_frame_fixed_pose_parity=True),indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("sources","candidates","aurora","brd","out"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    check(a.sources,a.candidates,a.aurora,a.brd,a.out)
