"""Cross score every dry and OLC-restored pool with identical query topology."""
import argparse
import itertools
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import time

import numpy as np
from scipy.spatial.distance import cdist

from restore import BENCH, write_csv
sys.path.insert(0,str(BENCH/"ligand_hydroxyl"))
from model import atoms, ligand_sites, topology

ARMS=("reported","dry","olc")


def components(lipid, block, query):
    """Continuous analytic Vina contributions for sugar vs C/O-only OLC (before table lookup and curl)."""
    l=atoms(block)
    _,donors=ligand_sites(query)
    sugar=[(i,a) for i,a in l.items() if a.kind!="HD"]
    serial=[i for i,a in sugar]
    q=np.array([a.xyz for i,a in sugar])
    qo=np.array([a.kind=="OA" for i,a in sugar])
    qh=np.array([bool(donors[i]) for i in serial])
    # Every sp3 carbon in GLC/BGC is covalently attached to oxygen: Vina C_P.
    mol,_,_=topology(query)
    if set(a.kind for i,a in sugar)!={"C","OA"} or any(
            all(b.GetAtomicNum()!=8 for b in a.GetNeighbors()) for a in mol.GetAtoms() if a.GetAtomicNum()==6):
        raise ValueError("The analytic decomposition is specific to these two sugar input graphs")
    out=[]
    for seq in (501,502,503):
        selected=[a for a in lipid.values() if int(a.residue[5:9])==seq]
        heavy=[a for a in selected if a.kind!="HD"]
        x=np.array([a.xyz for a in heavy])
        o=np.array([a.kind=="OA" for a in heavy])
        oh=np.array([min(np.linalg.norm(a.xyz-h.xyz) for h in selected if h.kind=="HD")<1.2 for a in heavy])
        r=cdist(q,x)
        gap=r-np.where(qo,1.7,1.9)[:,None]-np.where(o,1.7,1.9)[None,:]
        within=r<8
        hb=(qo[:,None]&oh[None,:])|(qh[:,None]&o[None,:])
        parts=dict(gauss1=-.035579*np.exp(-(gap/.5)**2),gauss2=-.005156*np.exp(-((gap-3.)/2.)**2),
                   repulsion=.840245*np.minimum(0.,gap)**2,
                   hydrophobic=np.zeros_like(r),hbond=-.587439*np.clip(-gap/.7,0,1)*hb)
        row={k:float((v*within).sum()) for k,v in parts.items()}
        ii,jj=np.unravel_index(np.argmin(r),r.shape)
        out.append(dict(seq=seq,nearest_sugar=serial[ii],nearest_lipid=heavy[jj].name,
                        distance_A=float(r.min()),within_cutoff=int(within.sum()),**row,
                        total=sum(row.values())))
    return out


def crystal_block(query,reference,serial):
    xyz,owners=ligand_sites(query)
    pos={i:np.array(v) for i,v in zip(serial,reference,strict=True)}
    for i,hs in owners.items():
        for h in hs:
            pos[h]=pos[i]+xyz[h].xyz-xyz[i].xyz
    lines=[]
    for l in query.splitlines():
        if l.startswith(("ATOM  ","HETATM")):
            v=pos[int(l[6:11])]
            l=l[:30]+"".join(f"{x:8.3f}" for x in v)+l[54:]
        lines.append(l)
    return "\n".join(lines)+"\n"


def cli_score(binary,rec,lig,pose,extra=()):
    cmd=[str(binary),"--receptor",str(rec),"--ligand",str(lig),"--pose_coordinates",str(pose),
         "--net_polar","native","--score_only","--cpu","1","--seed","43103","--verbosity","0",*extra]
    log=subprocess.check_output(cmd,text=True)
    native=next(l for l in log.splitlines() if l.startswith("Net polar native:")).split()[3:]
    return tuple(map(float,native))


def evaluate(w,o,binary):
    meta=json.loads((w/"metadata.json").read_text())
    lipid=atoms((o/"lipids.pdbqt").read_text())
    poses,cells,energy=[],[],[]
    start=time.perf_counter()
    with TemporaryDirectory() as tmp:
        pose=Path(tmp)/"frozen.pdbqt"
        for code in ("GLC","BGC"):
            m=meta["4ZW9-"+code]
            lig=w/"4ZW9"/(code+".pdbqt")
            query=lig.read_text()
            for pool,cap,seed in [("crystal",0,0),*itertools.product(("dry","olc"),(3000,12000),(43103,43117))]:
                if pool=="crystal":
                    blocks=[crystal_block(query,m["mapped"][0],m["serial"])]
                else:
                    suff="pro" if pool=="dry" else "olc"
                    path=(w if pool=="dry" else o)/"runs"/f"4ZW9-{code}_{cap}_{seed}_{suff}.pdbqt"
                    blocks=[b+"ENDMDL\n" for b in path.read_text().split("ENDMDL") if b.strip()]
                rr=[]
                for rank,b in enumerate(blocks,1):
                    common=dict(code=code,pool=pool,cap=cap,seed=seed,rank=rank)
                    coords=atoms(b)
                    d=np.asarray(m["mapped"])-np.array([coords[i].xyz for i in m["serial"]])
                    rmsd=float(np.sqrt((d*d).sum(axis=2).mean(axis=1)).min())
                    pose.write_text(b)
                    a=cli_score(binary,w/"4ZW9/receptor.pdbqt",lig,pose)
                    z=cli_score(binary,o/"4ZW9-OLC.pdbqt",lig,pose)
                    if a[1]!=z[1]:
                        raise ValueError("Intramolecular energy changed on same coordinates")
                    terms=components(lipid,b,query)
                    energy.extend(dict(**common,**t) for t in terms)
                    native=float(next(l for l in b.splitlines() if l.startswith("REMARK VINA RESULT:")).split()[3]) if pool!="crystal" else None
                    rr.append(dict(**common,rmsd=rmsd,reported=native,dry=a[2],olc=z[2],dry_inter=a[0],
                                   olc_inter=z[0],intra=a[1],raw_inter_delta=z[0]-a[0],
                                   continuous_lipid_total=sum(t["total"] for t in terms)))
                poses+=rr
                if pool=="crystal":
                    continue
                common={k:v for k,v in common.items() if k!="rank"}
                oracle=min(row["rmsd"] for row in rr)
                for arm in ARMS:
                    order=sorted(rr,key=lambda row:(row[arm],row["rank"]))
                    cells.append(dict(**common,score=arm,modes=len(rr),rank1=order[0]["rank"],
                        top1=order[0]["rmsd"],top5=min(t["rmsd"] for t in order[:5]),oracle=oracle))
    write_csv(o/"poses.csv",poses)
    write_csv(o/"cells.csv",cells)
    write_csv(o/"components.csv",energy)
    (o/"paired_cost.json").write_text(json.dumps(dict(wall_s=time.perf_counter()-start,cli_calls=2*len(poses)),indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("sugar",type=Path)
    p.add_argument("out",type=Path)
    p.add_argument("binary",type=Path)
    a=p.parse_args()
    evaluate(a.sugar.resolve(),a.out.resolve(),a.binary.resolve())
