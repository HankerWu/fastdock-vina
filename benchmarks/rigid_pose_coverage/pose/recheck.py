"""Independent coordinate replay, including feasible correspondences discarded by every gate."""
import argparse
import csv
from pathlib import Path
import sys
import tarfile

import numpy as np
from scipy.spatial import cKDTree

import run_panel as panel

sys.path[:0]=[str(panel.BENCH/"cavity_field_rank/audit"),str(panel.BENCH/"fragment_anchor/audit")]
import direct_space as independent
from edge_audit import align, raw_source, fits
from pocket_triangle import Pocket, _ligand_triangle
sys.path.insert(0,str(panel.HERE.parent/"engine"))
from coverage_gate import candidates


def check(p,out):
    components={"3ANR":"HRM","7AKL":"RK5","3FMZ":"2T1","1RBP":"RTL"}
    with tarfile.open(panel.BENCH/"cubic_grid/human_panel/depositions.tar.xz") as tar:
        for pdb,c in components.items():
            for file in (pdb+".pdb",c+".cif"): tar.extract("raw/"+file,p/"recheck",filter="data")
    frames={pdb:independent.reference(p/"recheck",pdb,c) for pdb,c in components.items()}
    independent.PDB.update({q:(query,rec) for q,(query,rec,_) in panel.CASES.items()})
    mode_rows=panel.old.csvrows(out/"modes.csv")
    selected={v["case"]:v["selected_stem"] for v in panel.old.csvrows(p/"extra_trajectory_selected.csv")}
    for row in panel.old.csvrows(out/"extra_trajectory_modes.csv"):
        name,seed=row["case"].split("_")
        mode_rows.append(dict(stem=selected[row["case"]],drug=name,seed=seed,arm="pro-extra",
                              rank=row["rank"],rmsd_A=row["rmsd_A"]))
    ranks={v["stem"]:v for v in panel.old.csvrows(out/"rank_stages.csv")}
    gates=[]
    for name,c in panel.old.configs(p).items():
        serial,target=independent.equivalent_targets(name,c["ligand"],frames)
        pocket=Pocket.from_pdbqt(c["receptor"],c["center"],c["size"])
        for seed in panel.CASES[name][2]:
            base=f"{name}_{seed}"
            ids,x=raw_source(p/"dump"/(base+"_edge")/"source.xyz")
            ref=target[:,[serial.index(s) for s in ids]]
            a=x[_ligand_triangle(x)]
            y=pocket.sites[pocket.triangles]
            rot=align(a,y)
            centers=((x.mean(0)-a.mean(0)) @ rot)+y.mean(1)
            pos=(x-a.mean(0))@rot + y.mean(1)[:,None,:]
            # Full atomic displacements are evaluated below for shortlisted nearest examples.
            lam,axis=np.linalg.eigh((x-x.mean(0)).T@(x-x.mean(0))/len(x))
            z=((rot.transpose(0,2,1)@axis)*np.sqrt(np.maximum(lam,0))).reshape(len(rot),9)
            z=np.column_stack((centers,z))
            ds=independent.distances(pos,ref)
            independent_feasible=((pocket.center-pocket.half_width<=pos)&(pos<=pocket.center+pocket.half_width)).all(axis=(1,2))
            eligible=candidates(pocket,x,ids)["eligible"]
            for arm in panel.ARMS:
                stem=base+"_"+arm
                with np.load(p/"dump"/stem/"bank.npz") as bank:
                    selected=bank["triangle_id"]
                    np.testing.assert_allclose(centers[selected],bank["center"],atol=1e-10,rtol=0)
                rr=ranks[stem]
                assert abs(ds[selected].min()-float(rr["generation_best_A"]))<1e-9
                distances,nn=cKDTree(z[selected]).query(z[eligible])
                pr=np.arange(0,len(distances),max(1,len(distances)//128))
                dim=np.linalg.norm((pos[eligible][pr]-pos[selected[nn[pr]]]).reshape(len(pr),-1),axis=1)/np.sqrt(len(x))
                np.testing.assert_allclose(distances[pr],dim,rtol=0,atol=1e-9)
                nearest=pocket.tree.query(pos[selected])[0].min(1)
                scores=panel.old.csvrows(p/"dump"/stem/"scores.csv")
                v=np.array([float(r["score"]) for r in scores])
                gates.append(dict(stem=stem,domain=len(ds),feasible=int(eligible.sum()),independent_feasible=int(independent_feasible.sum()),
                        support_disagreements=int((independent_feasible!=eligible).sum()),unfiltered_min_A=ds.min(),
                        unfiltered_le2=int((ds<=2).sum()),feasible_min_A=ds[eligible].min(),
                        feasible_le2=int((ds[eligible]<=2).sum()),rigid_lower_A=fits(x,ref,np.arange(len(x))),
                        gate_min_A=ds[selected].min(),gate_le2=int((ds[selected]<=2).sum()),
                        outside_box=int((~eligible[selected]).sum()),poses_with_contact_lt2A=int((nearest<2).sum()),
                        field_tied_minima=int((v==v.min()).sum()),
                        **dict(zip(("distance_median_A","distance_p90_A","distance_p95_A","distance_max_A"),
                                   np.percentile(distances,[50,90,95,100]),strict=True))))
            for mode in [m for m in mode_rows if m["drug"]==name and int(m["seed"])==seed]:
                stem=mode["stem"]
                if mode["arm"]=="pro-cost": stem+="-count"
                pool=(p/"raw"/(stem+".pdbqt")).read_text()
                block=pool.split("ENDMDL")[int(mode["rank"])-1]
                coords={int(l[6:11]):[float(l[j:j+8]) for j in (30,38,46)] for l in block.splitlines()
                        if l.startswith(("HETATM","ATOM  "))}
                error=np.linalg.norm(target-np.array([coords[s] for s in serial]),axis=(1,2)).min()/np.sqrt(len(serial))
                np.testing.assert_allclose(error,float(mode["rmsd_A"]),rtol=0,atol=1e-9)
    panel.table(out/"gates.csv",gates)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    for key in ("work","out"): cli.add_argument(key,type=lambda s:Path(s).resolve())
    a=cli.parse_args()
    check(a.work,a.out)
