"""Independent endpoint, chemistry, identity and raw AD4 map checks from written pools."""
import argparse
import csv
import itertools
import json
from pathlib import Path
import sys
import tarfile

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree
from scipy.stats import kendalltau, spearmanr

from assay_runner import ARMS, LIVE, RECEPTORS, ROOT, SEED, SERIES, prep, sha, table
sys.path.insert(0, str(ROOT/"benchmarks/cavity_field_rank/fkbp"))
from frozen_assay import old
from check_assay import rms


def trace(stem, native, serials, blocks, key):
    h, line=stem.with_suffix(".source").read_text().splitlines()
    raw_order=list(map(int,h[1:].split()))
    reorder=np.argsort(raw_order)
    source=np.array(line.split()[3:],float).reshape(-1,3)[reorder]
    sorted_serial=np.array(raw_order)[reorder].tolist()
    chemical=[sorted_serial.index(i) for i in serials]
    snap_order=[raw_order.index(i) for i in serials]
    ranklog=json.loads(stem.with_suffix(".rank").read_text())
    with np.load(stem.with_suffix(".npz")) as bank:
        poses=(source-source.mean(0)) @ bank["rotation"].transpose(0,2,1)+bank["center"][:,None]
        idx=np.argsort(bank["true_score"],kind="stable")
        i=int(idx[0])
        assert i==ranklog["winner"] and len(idx)==1024
        scores=bank["true_score"]
        d=rms(poses[:,chemical],native) if native is not None else None
        anchor=" ".join(map(str,bank["anchor_serials"]))
        component=" ".join(map(str,bank["component_serials"]))
        fallback=int(bank["fallback"])
    snapshots=[]
    for l in stem.with_suffix(".xyz").read_text().splitlines()[1:]:
        name,tr,pos,text=l.split(" ",3)
        snapshots.append((name,int(tr),int(pos),np.fromstring(text,sep=" ").reshape(-1,3)[snap_order]))
    history={int(l.split()[0]):list(map(int,l.split()[4:])) for l in stem.with_suffix(".history").read_text().splitlines()}
    stages={k:x for k,tr,pos,x in snapshots if k in ("selected","hunt","second")}
    np.testing.assert_allclose(stages["selected"],poses[i,chemical],rtol=0,atol=1e-9)
    at_end=[p for k,tr,p,x in snapshots if k=="end_trajectory" and tr==1 and np.array_equal(x,stages["second"])]
    merged=[p for k,tr,p,x in snapshots if k=="merged_ad4" and np.array_equal(x,stages["second"])]
    final=[j for j,b in enumerate(blocks,1) if np.max(abs(old.xyz(b,serials)-stages["second"]))<=.00051]
    assert bool(1 in history[max(history)])==bool(at_end) and (not final or merged) and (not merged or at_end)
    regens=[p for k,tr,p,x in snapshots if native is not None and k=="offered" and tr==1 and p>1 and rms(x,native)<=2]
    near=(idx+1) if d is None else np.flatnonzero(d[idx]<=2)+1
    return dict(**key, index=i,score=scores[i],anchor_serials=anchor,component_serials=component,fallback=fallback,
        selected_A=d[i] if d is not None else "",best_generated_A=min(d) if d is not None else "",
        generated_le2=int(sum(d<=2)) if d is not None else "",first_le2_rank=near[0] if d is not None and len(near) else "",
        hunt_A=float(rms(stages["hunt"],native)) if native is not None else "",
        second_A=float(rms(stages["second"],native)) if native is not None else "",entered=int(1 in history[1]),
        evicted_step=next((k for k,v in history.items() if 1 not in v),""),
        end_rank=next(iter(at_end),""),merge_rank=next(iter(merged),""),released_rank=next(iter(final),""),
        first_later_offered_le2=next(iter(regens),""), **ranklog)


def check(work,out):
    out.mkdir(parents=True,exist_ok=True)
    input_hashes=json.loads((work/"inputs.json").read_text())
    for p,d in input_hashes["files"].items(): assert sha(work/p)==d
    dep=SERIES/"human_panel/depositions.tar.xz"
    assay=SERIES/"human_panel/affinity/ligands.csv"
    with tarfile.open(dep) as tar:
        for f in ("1FKG.pdb","1FKH.pdb","SB3.cif","SBX.cif"): tar.extract("raw/"+f,work,filter="data")
    refs={r:old.reference(work,r,c) for r,c in zip(RECEPTORS,("SB3","SBX"),strict=True)}
    with assay.open() as f: labels={int(r["compound"]):r for r in csv.DictReader(f)}
    # Direct Table II transcription, JACS 115 (1993) p. 9927.
    measured=[(186,186),(110,110),(12,12),(250,250),(10,10),(7,7),(300,600),(300,300)]
    assert measured==[(float(labels[i]["ki_low_nm"]),float(labels[i]["ki_high_nm"])) for i in range(4,12)]
    (out/"inputs.json").write_text(json.dumps(dict(**input_hashes,evaluation={str(p.relative_to(ROOT)):sha(p) for p in (dep,assay)}),indent=2)+"\n")
    runs,cells,modes,traces=[],[],[],[]
    gate_rows=[]
    top,rebuilt,desc={},{},{}
    for rec in RECEPTORS:
        tree=cKDTree([[float(l[k:k+8]) for k in (30,38,46)] for l in (work/rec/"receptor.pdbqt").read_text().splitlines()
                      if l.startswith(("ATOM","HETATM")) and l.split()[-1] not in ("H","HD","HS")])
        fields={s:old.read_map(work/"maps"/rec/f"receptor.{s}.map") for s in ("A","C","N","OA","e","d")}
        for n in labels:
            lig=(work/"affinity"/f"{n}.pdbqt").read_text()
            mol,serials,bonds=old.molecule(lig)
            canonical=Chem.MolToSmiles(mol)
            assert canonical==Chem.MolToSmiles(Chem.MolFromSmiles(labels[n]["smiles"]))
            fn=lambda x:np.linalg.norm(x[bonds[:,0]]-x[bonds[:,1]],axis=-1)
            lengths=fn(old.xyz(lig,serials))
            crystal={8:"1FKG",9:"1FKH"}.get(n)
            native=old.crystal_positions(refs[crystal],refs[rec],mol) if crystal else None
            desc[n]=dict(compound=n,heavy_atoms=len(serials),charge=Chem.GetFormalCharge(mol),torsdof=int(lig.split("TORSDOF ")[-1].split()[0]))
            for arm in ARMS:
                key=dict(arm=arm,receiver=rec,compound=n)
                stem=work/"pools"/arm/f"{rec}_{n}"
                row=json.loads(stem.with_suffix(".timing").read_text())
                assert row["sha256"]==sha(stem.with_suffix(".pdbqt")) and row["seed"]==SEED
                assert f"(random seed: {SEED})" in stem.with_suffix(".log").read_text()
                blocks=[b for b in stem.with_suffix(".pdbqt").read_text().split("ENDMDL") if b.strip()]
                assert len(blocks)==row["modes"]
                runs.append(row)
                errors=[]
                for rank,b in enumerate(blocks,1):
                    x=old.xyz(b,serials)
                    m,s,e=old.molecule(b)
                    assert s==serials and np.array_equal(e,bonds) and Chem.MolToSmiles(m)==canonical
                    pred=float(next(l.split()[3] for l in b.splitlines() if l.startswith("REMARK VINA RESULT:")))
                    recalc=sum(old.calculate_components(b.splitlines(),fields).values())
                    err=float(rms(x,native)) if crystal else ""
                    if crystal: errors.append(err)
                    modes.append(dict(**key,rank=rank,score=pred,from_xyz=recalc,rmsd_A=err,min_contact_A=tree.query(x)[0].min(),
                                      stereo_changed=int(old.stereo(m,x)!=canonical),max_bond_delta_A=np.max(np.abs(fn(x)-lengths))))
                    if rank==1: top[arm,rec,n],rebuilt[arm,rec,n]=pred,recalc
                if crystal and errors:
                    cells.append(dict(**key,top1=errors[0],best5=min(errors[:5]),oracle=min(errors),regret=errors[0]-min(errors)))
                if arm in LIVE: traces.append(trace(stem,native,serials,blocks,key))
            if crystal:
                sys.path[:0]=[str(Path(__file__).resolve().parents[1]/"engine"),str(ROOT/"tools")]
                from anchors import all_transforms
                from pocket_triangle import Pocket
                with (work/"maps"/rec/"receptor.C.map").open() as h:
                    headers=[next(h).split() for _ in range(6)]
                pocket=Pocket.from_pdbqt(work/rec/"receptor.pdbqt", np.array(headers[5][1:],float),
                                         np.array(headers[4][1:],int)*float(headers[3][1]))
                for arm in ("whole","fragment"):
                    stem=work/"pools"/arm/f"{rec}_{n}"
                    h, line=stem.with_suffix(".source").read_text().splitlines()
                    raw_order=np.array(h[1:].split(),int)
                    order=np.argsort(raw_order)
                    src=np.array(line.split()[3:],float).reshape(-1,3)[order]
                    chem=[raw_order[order].tolist().index(i) for i in serials]
                    with np.load(stem.with_suffix(".npz")) as b:
                        all_fits=all_transforms(pocket,src,b["anchor"])
                        ranked=np.argsort(np.round(all_fits["edge_residual_A2"],10),kind="stable")
                        assert np.array_equal(all_fits["site_triangles"][ranked[:1024]],b["site_triangles"])
                    errs=[]
                    for chunk in np.array_split(ranked,128):
                        fit=(src-src.mean(0)) @ all_fits["matched_rotation"][chunk].transpose(0,2,1)+all_fits["center"][chunk,None]
                        errs.extend(rms(fit[:,chem],native))
                    errs=np.asarray(errs)
                    pres=np.flatnonzero(errs<=2)
                    gate_rows.append(dict(receiver=rec,compound=n,arm=arm,all=len(errs),prefilter_min_A=float(errs.min()),
                        postfilter_min_A=float(errs[:1024].min()),first_le2_edge_rank=int(pres[0]+1) if len(pres) else "",
                        pre_le2=len(pres),post_le2=int((errs[:1024]<=2).sum())))
            for ext in (".source",".q"):
                assert len({sha((work/"pools"/a/f"{rec}_{n}").with_suffix(ext)) for a in LIVE})==1
    assayrows,pairs,scores,quantiles,costs=[],[],[],[],[]
    for a in ARMS:
        for n in labels:
            for d in (top,rebuilt): d[a,"mean",n]=float(np.mean([d[a,r,n] for r in RECEPTORS]))
        for rec in (*RECEPTORS,"mean"):
            point=[i for i in labels if measured[i-4][0]==measured[i-4][1]]
            pred=[top[a,rec,i] for i in point]
            obs=[measured[i-4][0] for i in point]
            credits=[]
            for i,j in itertools.permutations(labels,2):
                fold=measured[j-4][0]/measured[i-4][1]
                if fold<3: continue
                margin=top[a,rec,j]-top[a,rec,i]
                xyz_margin=rebuilt[a,rec,j]-rebuilt[a,rec,i]
                credit=(np.sign(margin)+1)/2
                credits.append(credit)
                pairs.append(dict(arm=a,receiver=rec,stronger=i,weaker=j,minimum_fold=fold,weak_minus_strong=margin,
                                  xyz_margin=xyz_margin,credit=credit,xyz_credit=(np.sign(xyz_margin)+1)/2))
            assert len(credits)==15
            scores.extend(dict(arm=a,receiver=rec,compound=n,score=top[a,rec,n],from_xyz=rebuilt[a,rec,n],ki_low=measured[n-4][0],
                               ki_high=measured[n-4][1]) for n in labels)
            assayrows.append(dict(arm=a,receiver=rec,correct=sum(credits),denominator=len(credits),points=len(point),
                rho=spearmanr(obs,pred).statistic,tau_b=kendalltau(obs,pred).statistic,
                cliff8_10=top[a,rec,10]-top[a,rec,8],cliff9_11=top[a,rec,11]-top[a,rec,9]))
        for n in ("all",8,9):
            rows=[r for r in cells if r["arm"]==a and (n=="all" or r["compound"]==n)]
            x=np.array([r["top1"] for r in rows])
            quantiles.append(dict(arm=a,compound=n,cells=len(x),mean=x.mean(),**dict(zip(("q25","median","q75","p90","p95"),
                             np.quantile(x,[.25,.5,.75,.9,.95]),strict=True)),**{k:np.mean([r[k] for r in rows]) for k in ("best5","oracle","regret")},
                             gt5=sum(x>5),gt10=sum(x>10)))
        rows=[r for r in runs if r["arm"]==a]
        m=[r for r in modes if r["arm"]==a]
        costs.append(dict(arm=a,cells=len(rows),empty=sum(r["modes"]==0 for r in rows),modes=len(m),
            stereo_changes=sum(r["stereo_changed"] for r in m),max_bond_delta_A=max(r["max_bond_delta_A"] for r in m),
            min_heavy_contact_A=min(r["min_contact_A"] for r in m),max_score_xyz_delta=max(abs(r["score"]-r["from_xyz"]) for r in m),
            **{k:sum(r[k] for r in rows) for k in ("hunt_evals","second_evals","coarse_channels","downstream_channels",
                                                  "wall_s","cpu_s","official_counter_replay_wall_s","official_counter_replay_cpu_s")}))
    with (work/"unmodified_pro.csv").open() as f: production=list(csv.DictReader(f))
    for r in production:
        assert r["sha256"]==sha(work/"pools"/r["arm"]/(r["receiver"]+"_"+r["compound"]+".pdbqt"))
    biases=[]
    for b in ("base","whole","sham"):
        delta=[top["fragment","mean",n]-top[b,"mean",n] for n in labels]
        biases.append(dict(arm="fragment",control=b,n=len(delta),charge_levels=len({r["charge"] for r in desc.values()}),
                           **{"rho_"+k:spearmanr(delta,[desc[n][k] for n in labels]).statistic for k in ("heavy_atoms","torsdof")}))
    for name,data in (("assay",assayrows),("scores",scores),("pairs",pairs),("runs",runs),("pose_cells",cells),
                      ("modes",modes),("stages",traces),("distributions",quantiles),("cost",costs),("bias",biases),
                      ("descriptors",list(desc.values())),("gate",gate_rows),("unmodified_pro",production)):
        table(out/(name+".csv"),data)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    check(a.work.resolve(),a.out.resolve())
