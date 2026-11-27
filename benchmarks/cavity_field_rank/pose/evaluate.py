"""Evaluate released pools and single-event lineage in the receiver frame after predictions are locked."""
import argparse
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import sys
import tarfile

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

SPEC = importlib.util.spec_from_file_location("field_panel", Path(__file__).with_name("run.py"))
panel = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(panel)
from pose_metrics import inspect_fork, scan
from evaluate_acid import xyz, neutral_mol, reference, transform
from pose_runs import table, digest

COMPONENTS = {"1RBP":"RTL","3FMZ":"2T1","4GHI":"0X3","5TBM":"79A"}

def reference_inputs(p):
    files = {key:sorted((p/"inputs").glob(key+"-ligand.pdbqt"))[0] for key in panel.CASES}
    data = {}
    with tarfile.open(panel.BENCH/"cubic_grid/human_panel/depositions.tar.xz") as f:
        for pdb,component in COMPONENTS.items():
            for element in (pdb+".pdb",component+".cif"):
                f.extract("raw/"+element,p/"refs",filter="data")
    refs={key:reference(p/"refs",key,component) for key,component in COMPONENTS.items()}
    # Primary BRD9 correspondence from SIFTS and CCD in the deposited fixture; audit can refit directly to CA.
    with tarfile.open(panel.LIB["brd"]) as f:
        meta=json.load(f.extractfile("metadata.json"))
    cross=scan.load_path("cross",panel.BENCH/"paired_crossdocking/run.py")
    torsion=scan.load_path("torsions",panel.BENCH/"torsion_proposals/paired.py")
    for name,(q,r,_) in panel.CASES.items():
        lines=files[name].read_text().splitlines()
        mol=Chem.MolFromSmiles(next(l.split()[2] for l in lines if l.startswith("REMARK SMILES ")))
        mapping={}
        for l in lines:
            if l.startswith("REMARK SMILES IDX "):
                v=list(map(int,l.split()[3:]))
                mapping.update(zip(v[::2],v[1::2],strict=True))
        ids=[mapping[i+1] for i in range(mol.GetNumAtoms())]
        if name=="82I":
            a,b=meta[q],meta[r]
            rot,t,_,_=cross.alignment(*torsion.aligned_receivers(a,b))
            exp=np.asarray(a["mapped"]) @ rot+t
            assert ids==a["serial"]
        else:
            ref=refs[q][2]
            site=transform(q,r,refs)(ref.GetConformer().GetPositions())
            neutral=neutral_mol(mol)
            maps=ref.GetSubstructMatches(neutral,uniquify=False,useChirality=True,maxMatches=10001)
            assert 0<len(maps)<10001
            exp=site[np.asarray(maps)]
            if name=="3FMZ-A":
                ((_, a,b),)=neutral.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]"))
                swap=exp.copy()
                swap[:,[a,b]]=exp[:,[b,a]]
                exp=np.concatenate([exp,swap])
        rec=(p/"inputs"/(r+".pdbqt")).read_text()
        heavy=[int(l[6:11]) for l in rec.splitlines() if l.startswith(("ATOM  ","HETATM")) and l.split()[-1] not in ("H","HD","HS")]
        data[name]=(dict(refs=exp,serials=ids), mol, xyz("\n".join(lines),ids), cKDTree(xyz(rec,heavy)))
    return data

def errors(x, refs):
    sq=np.square(x[...,None,:,:]-refs).sum(axis=-1).mean(axis=-1)
    return np.sqrt(sq.min(axis=-1))

def quality(x,src,mol,tree):
    e=np.array([(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in mol.GetBonds()])
    lens=lambda a:np.linalg.norm(a[e[:,0]]-a[e[:,1]],axis=1)
    inv=0
    stereo=Chem.FindMolChiralCenters(mol,includeUnassigned=False)
    for i,_ in stereo:
        neighbors=[v.GetIdx() for v in mol.GetAtomWithIdx(i).GetNeighbors()][:3]
        signed=lambda a:np.linalg.det(a[neighbors]-a[i])
        inv+=int(signed(x)*signed(src)<0)
    return dict(max_bond_drift_A=max(abs(lens(x)-lens(src))),defined_centres=len(stereo),
                inversions=inv,min_receptor_distance_A=tree.query(x)[0].min())

def analyze(p,out):
    out.mkdir(exist_ok=True,parents=True)
    geom=reference_inputs(p)
    costs={r["stem"]:r for n in ("control_cost","dock_cost") for r in panel.csvrows(p/(n+".csv"))}
    cells,modes,chem,stages,residency,forks,candidates=[],[],[],[],[],[],[]
    for name,(query,rec,seeds) in panel.CASES.items():
        cell,mol,src,tree=geom[name]
        group="BRD9" if name=="82I" else ("RBP4" if name.startswith("3FMZ") else "EPAS1")
        for seed in seeds:
            case=f"{name}_{seed}"
            bank=p/"banks"/case
            s=np.load(bank.with_suffix(".source.npz"))
            c=np.load(bank.with_suffix(".npz"))
            order=[list(s["serials"]).index(v) for v in cell["serials"]]
            x=s["xyz"]
            pos=(x-x.mean(0)) @ c["matched_rotation"].transpose(0,2,1)+c["center"][:,None,:]
            vals=errors(pos[:,order],cell["refs"])
            stats={}
            for arm in ("shape","true","perm","pro","pro-high","official","official-high"):
                stem=f"{case}_{arm}"
                key=dict(stem=stem,panel=group,query=query,receiver=rec,drug=name,seed=seed,arm=arm)
                raw=p/"raw"/stem
                pool=raw.with_suffix(".pdbqt")
                assert digest(pool)==costs[stem]["pool_sha256"]
                errs=[]
                for rank,block in enumerate(pool.read_text().split("ENDMDL")[:-1],1):
                    xyz_a=xyz(block,cell["serials"])
                    r=float(errors(xyz_a,cell["refs"]))
                    v=float(next(l.split()[3] for l in block.splitlines() if l.startswith("REMARK VINA RESULT:")))
                    errs.append(r)
                    modes.append(dict(**key,rank=rank,rmsd_A=r,score=v))
                    chem.append(dict(**key,rank=rank,**quality(xyz_a,src,mol,tree)))
                cells.append(dict(**key,n=len(errs),top1_A=errs[0] if errs else "",
                         best5_A=min(errs[:5],default=""),oracle_A=min(errs,default=""),
                         regret_A=errs[0]-min(errs) if errs else ""))
                if arm not in ("shape","true","perm"):
                    continue
                if arm=="shape":
                    scores=c["matched_score"]
                    ids=np.argsort(np.round(scores,10),kind="stable")
                    i=int(c["matched_best"])
                    nchanged,nquery="",""
                else:
                    rows=panel.csvrows(raw.with_suffix(".csv"))
                    assert [int(v["index"]) for v in rows]==list(range(1024))
                    scores=np.array([float(v["score"]) for v in rows])
                    ids=np.argsort(scores,kind="stable")
                    l=next(l for l in raw.with_suffix(".inject").read_text().splitlines() if l.startswith("field_rank "))
                    i=int(l.split()[3])
                    nchanged,nquery=(int(l.split()[k]) for k in (9,11))
                    assert nquery==len(x)*1024
                assert ids[0]==i
                ranks=np.argsort(ids)+1
                stats[arm]=(scores,ranks)
                st,hist,chain=inspect_fork(p,stem,cell,1,1)
                stages+=st
                residency+=hist
                true_order=np.array(list(map(int,raw.with_suffix(".stages").read_text().splitlines()[0][1:].split())))
                np.testing.assert_array_equal(true_order,s["serials"])
                selected_xyz=np.array(next(l.split()[3:] for l in raw.with_suffix(".stages").read_text().splitlines() if l.startswith("selected ")),float).reshape(-1,3)
                np.testing.assert_allclose(selected_xyz,pos[i],atol=1e-9,rtol=0)
                second=next((v["rmsd_A"] for v in st if v["stage"]=="second"),"")
                forks.append(dict(**key,selected_index=i,heavy_atoms=len(x),n_type_changed=nchanged,rank_field_queries=nquery,
                    n_le2=int(sum(vals<=2)),n_le3=int(sum(vals<=3)),generation_best_A=vals.min(),
                    first_le2_rank=min(ranks[vals<=2],default=""),first_le3_rank=min(ranks[vals<=3],default=""),
                    shape_rank_of_selected=int(stats["shape"][1][i]),source_A=float(errors(x[order],cell["refs"])),second_A=second,
                    first_later_offered_le2_step=next((v["index"] for v in st if v["stage"]=="offered" and v["scope"]==1 and v["index"]>1 and v["rmsd_A"]<=2),""),
                    **{k:v for k,v in chain.items() if k!="stem"}))
            for i,v in enumerate(vals):
                row=dict(case=case,index=i,rmsd_A=v)
                for arm,(scores,ranks) in stats.items():
                    row[arm+"_score"]=scores[i]
                    row[arm+"_rank"]=ranks[i]
                candidates.append(row)
    for name,rows in (("cells",cells),("modes",modes),("chemistry",chem),("stages",stages),("residency",residency),("rank_stages",forks)):
        table(out/(name+".csv"),rows)
    with gzip.open(out/"candidates.csv.gz","wt",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=list(candidates[0]));w.writeheader();w.writerows(candidates)
    dist=[]
    for arm in sorted({c["arm"] for c in cells}):
        for group in ("all","BRD9","EPAS1","RBP4",*geom):
            subset=[r for r in cells if r["arm"]==arm and (group=="all" or group in (r["drug"],r["panel"]))]
            vals=[r for r in subset if r["n"]]
            top=np.array([r["top1_A"] for r in vals])
            dist.append(dict(group=group,arm=arm,cells=len(subset),empty=len(subset)-len(vals),
                **{k:np.mean([r[k] for r in vals]) for k in ("top1_A","best5_A","oracle_A","regret_A")},
                **dict(zip(("q25","median","q75","p90","p95"),np.percentile(top,[25,50,75,90,95]),strict=True)),
                gt5=sum(top>5),gt10=sum(top>10)))
    table(out/"distributions.csv",dist)
    panel.write_json(out/"references.json",dict(
        depositions={str(f.relative_to(panel.ROOT)):digest(f) for f in (panel.LIB["brd"],panel.BENCH/"cubic_grid/human_panel/depositions.tar.xz")},
        atom_serials_and_equivalent_mappings={k:dict(serials=v[0]["serials"],maps=len(v[0]["refs"])) for k,v in geom.items()}))

if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("work",type=Path); a.add_argument("output",type=Path)
    args=a.parse_args()
    analyze(args.work.resolve(),args.output.resolve())
