"""Post-outcome label-only localization of losses at the 1,024-triangle gate."""

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist

from evaluate_candidates import study
from generate import Pocket, write


def run(w, dest):
    systems = study.scan.config(w / "aurora", w / "brd", dest)
    summaries, costs, pockets, raw = [], [], {}, {}
    for row in csv.DictReader((w / "sources.csv").open()):
        source, drug = row["source"], row["drug"]
        with np.load(row["npz"]) as parent, np.load(w / "candidates" / (source + ".npz")) as observed:
            x = parent["xyz"]
            serials = parent["serials"].tolist()
            atoms = observed["ligand_triangle"]
            keep = observed["site_triangles"]
        k = row["receptor"], tuple(float(row["center_" + i]) for i in "xyz"), float(row["size"])
        if k not in pockets:
            pockets[k] = Pocket.from_pdbqt(*k)
        pocket = pockets[k]
        t0, c0 = time.perf_counter(), time.process_time()
        refs = systems[drug]["refs"][:, [systems[drug]["serials"].index(i) for i in serials]]
        src = x[atoms]
        priority = np.argsort(np.round(((pocket.edges - pdist(src))**2).sum(axis=1),10),kind="stable")
        assert np.array_equal(pocket.triangles[priority[:len(keep)]],keep)
        nearest = []
        for a in range(0,len(priority),256):
            triangles = pocket.triangles[priority[a:a+256]]
            dest_xyz = pocket.sites[triangles]
            mean = dest_xyz.mean(axis=1)
            centered = dest_xyz - mean[:,None,:]
            u,_,vt = np.linalg.svd((src-src.mean(axis=0)).T @ centered)
            u[:,:,-1] *= np.linalg.det(u @ vt)[:,None]
            poses = (x-src.mean(axis=0)) @ (u @ vt) + mean[:,None,:]
            delta = poses[:,None,:,:]-refs[None,:,:,:]
            nearest.extend(np.sqrt(np.square(delta).sum(axis=-1).mean(axis=-1).min(axis=1)))
        nearest = np.asarray(nearest)
        raw[source] = nearest.astype(np.float32)
        best = nearest.argmin()
        r = dict(source=source, drug=drug, seed=row["seed"], population=len(priority),
                 **{f"site{i+1}":j for i,j in enumerate(pocket.triangles[priority[best]])},
                 before_gate_min_A=nearest[best], best_edge_rank=best+1,
                 after_gate_min_A=nearest[:1024].min(), kept_le2=(nearest[:1024]<=2).sum())
        for t in (2,3):
            ids = np.flatnonzero(nearest<=t)
            r[f"all_le{t}"] = len(ids)
            r[f"first_edge_rank_le{t}"] = ids[0]+1 if len(ids) else None
        summaries.append(r)
        costs.append(dict(source=source,wall_s=time.perf_counter()-t0,cpu_s=time.process_time()-c0))
    write(dest/"posthoc_edge_gate.csv",summaries)
    write(dest/"posthoc_cost.csv",costs)
    np.savez_compressed(dest/"posthoc_rmsd_by_edge_rank.npz",**raw)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    run(a.work,a.out)
