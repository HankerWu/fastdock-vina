"""Predeclared label-assisted forks. Freeze the 48 blind pools before use."""
import argparse
from pathlib import Path
import json
import sys

import numpy as np
from scipy.spatial.transform import Rotation

from edge_audit import frames, raw_source, equivalent_targets, table, TRIPLET
sys.path.insert(0,str(TRIPLET))
from pose_runs import invoke
from pose_metrics import inspect_fork

CASES=('3FMZ-A_93241','3FMZ-HA_93241','4GHI_93257','82I_91813','82I_91823')


def read_pool(path,ids,ref):
    errors=[]
    for block in path.read_text().split('ENDMDL')[:-1]:
        coords={int(l[6:11]):[float(l[k:k+8]) for k in (30,38,46)]
                for l in block.splitlines() if l.startswith(('ATOM  ','HETATM'))}
        x=np.array([coords[s] for s in ids])
        errors.append(float(np.linalg.norm(ref-x,axis=(1,2)).min()/np.sqrt(len(x))))
    return errors


def run(p,w,observer,out):
    if sum(len(list((p/'raw').glob('*_'+arm+'.pdbqt'))) for arm in ('whole','fragment','sham'))!=48:
        raise ValueError('incomplete blind predictions')
    (w/'raw').mkdir(parents=True,exist_ok=True)
    geo=frames(w)
    configs=json.loads((p/'inputs.json').read_text())['cells']
    rows,choices,costs=[],[],[]
    for case in CASES:
        name,seed=case.split('_')
        cfg=configs[name].copy()
        for k in ('ligand','bonds','receptor'):
            cfg[k]=(p/cfg[k]).resolve()
        order,ref=equivalent_targets(name,cfg['ligand'],geo)
        for arm in ('whole','fragment'):
            stem=f'{case}_{arm}'
            src=p/'dump'/stem
            ids,x=raw_source(src/'source.xyz')
            target=ref[:,[order.index(s) for s in ids]]
            q=np.loadtxt(src/'source.q')[2:]
            with np.load(src/'bank.npz') as bank:
                rot,centers=bank['matched_rotation'],bank['center']
            xyz=(x-x.mean(0))@rot.transpose(0,2,1)+centers[:,None,:]
            e=np.sqrt(np.min([np.square(xyz-r).sum(2).mean(1) for r in target],axis=0))
            scores=np.loadtxt(src/'scores.csv',delimiter=',',skiprows=1)[:,1]
            ranked=np.argsort(scores,kind='stable')
            choice=next(((j+1,i) for j,i in enumerate(ranked) if e[i]<=2.),None)
            choices.append(dict(stem=stem,eligible=int(sum(e<=2)),tested=int(choice is not None)))
            if choice is None:
                continue
            rank,i=choice
            rigid=Rotation.from_matrix(rot[i])
            pos=rigid.apply(q[:3]-x.mean(0))+centers[i]
            quat=(rigid*Rotation.from_quat(q[[4,5,6,3]])).as_quat()[[3,0,1,2]]
            base=w/'raw'/stem
            event=base.with_suffix('.input')
            event.write_text('1 1 '+' '.join(format(z,'.17g') for z in np.r_[pos,quat])+'\n')
            result=invoke(w,costs,stem,observer,cfg,int(seed),env=dict(FD_RIGID_INTERVENTION=event.resolve(),
                FD_RIGID_LOG=base.with_suffix('.inject').resolve(),FD_RIGID_COORDS=base.with_suffix('.stages').resolve(),
                FD_RIGID_HISTORY=base.with_suffix('.history').resolve()))
            stages,_,lineage=inspect_fork(w,stem,dict(refs=ref,serials=order),1,1)
            np.testing.assert_allclose(e[i],lineage['selected_A'],atol=1e-9)
            pool=read_pool(result,order,ref)
            rows.append(dict(case=case,arm=arm,candidate_index=i,field_rank=rank,
                second_A=next((s['rmsd_A'] for s in stages if s['stage']=='second'),''),
                top1_A=pool[0],best5_A=min(pool[:5]),oracle_A=min(pool),**lineage))
    table(out/'fork_eligibility.csv',choices)
    table(out/'near_forks.csv',rows)
    table(out/'fork_cost.csv',costs)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ('pose','work','observer','output'):
        parser.add_argument(k,type=lambda x:Path(x).resolve())
    a=parser.parse_args()
    run(a.pose,a.work,a.observer,a.output)
