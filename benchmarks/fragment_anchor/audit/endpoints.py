"""Independent receiver alignment, atom matching and recount from original outputs and Table II."""
import argparse
import csv
import gzip
import itertools
import json
from pathlib import Path
import time

import numpy as np
from scipy.stats import kendalltau, spearmanr

from edge_audit import frames, raw_source, equivalent_targets, table
from forks import read_pool

ARMS=('whole','fragment','sham','pro','pro-high','official','official-high')


def rows(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:
        return list(csv.DictReader(f))


def pose(p, work, results, out):
    frame=frames(work)
    config=json.loads((p/'inputs.json').read_text())['cells']
    cells={r['stem']:r for r in rows(results/'cells.csv')}
    modes={(r['stem'],int(r['rank'])):float(r['rmsd_A']) for r in rows(results/'modes.csv')}
    candidates={(r['stem'],int(r['index'])):float(r['rmsd_A']) for r in rows(results/'candidates.csv.gz')}
    stages={(r['stem'],r['stage'],r['scope'],r['index']):float(r['rmsd_A']) for r in rows(results/'stages.csv')}
    nm,nc,ns=len(modes),len(candidates),len(stages)
    evidence=[]
    err=dict(candidate_A=0.,stage_A=0.,release_A=0.,endpoint_A=0.)
    for stem,row in sorted(cells.items()):
        drug,seed,arm=stem.split('_')
        assert arm in ARMS
        ids,ref=equivalent_targets(drug,p/config[drug]['ligand'],frame)
        actual=read_pool(p/'raw'/(stem+'.pdbqt'),ids,ref)
        assert len(actual)==int(row['n'])
        for k,d in enumerate(actual,1):
            err['release_A']=max(err['release_A'],abs(modes.pop((stem,k))-d))
        if actual:
            e=[actual[0],min(actual[:5]),min(actual),actual[0]-min(actual)]
            err['endpoint_A']=max(err['endpoint_A'],max(abs(a-float(row[k])) for a,k in
                zip(e,('top1_A','best5_A','oracle_A','regret_A'),strict=True)))
        evidence.append(dict(stem=stem,n=len(actual),top1_A=actual[0] if actual else '',oracle_A=min(actual,default='')))
        if arm not in ('whole','fragment','sham'):
            continue
        d=p/'dump'/stem
        serial,x=raw_source(d/'source.xyz')
        reordered=ref[:,[ids.index(s) for s in serial]]
        with np.load(d/'bank.npz') as bank:
            xyz=(x-x.mean(0))@bank['matched_rotation'].transpose(0,2,1)+bank['center'][:,None,:]
        ds=np.min([np.linalg.norm(xyz-r,axis=(1,2)) for r in reordered],axis=0)/np.sqrt(len(x))
        err['candidate_A']=max(err['candidate_A'],max(abs(y-candidates.pop((stem,i))) for i,y in enumerate(ds)))
        with (d/'stages.xyz').open() as trace:
            order=list(map(int,next(trace)[1:].split()))
            targets=ref[:,[ids.index(s) for s in order]]
            for line in trace:
                tag,scope,index,coords=line.split(maxsplit=3)
                xyz=np.fromstring(coords,sep=' ').reshape(-1,3)
                e=np.min(np.linalg.norm(targets-xyz,axis=(1,2)))/np.sqrt(len(xyz))
                err['stage_A']=max(err['stage_A'],abs(e-stages.pop((stem,tag,scope,index))))
    assert not modes and not candidates and not stages and len(cells)==112
    assert max(err.values())<1e-8
    table(out/'released.csv',evidence)
    return dict(pools=len(cells),modes=nm,candidates=nc,stages=ns,max_difference=err)


def assay(p,out,submitted):
    # Independent Table II transcription: Holt et al., JACS 115, p. 9927. nM.
    ranges=[(186,186),(110,110),(12,12),(250,250),(10,10),(7,7),(300,600),(300,300)]
    reliable=[(i,j) for i,j in itertools.permutations(range(8),2) if ranges[j][0]>=3*ranges[i][1]]
    exact=[i for i,(a,b) in enumerate(ranges) if a==b]
    assert len(reliable)==15 and len(exact)==7
    true=np.array([ranges[i][0] for i in exact])
    records=[]
    for arm in ('official','base','whole','fragment','sham','high','official_high'):
        separate=[]
        for receptor in ('1FKG','1FKH'):
            score=[]
            for n in range(4,12):
                file=p/'pools'/arm/f'{receptor}_{n}.pdbqt'
                with file.open() as f:
                    score.append(float(next(line for line in f if line.startswith('REMARK VINA RESULT:')).split()[3]))
            separate.append(np.array(score))
        for rec,pred in zip(('1FKG','1FKH','mean'),[*separate,np.mean(separate,axis=0)],strict=True):
            wins=sum(1. if pred[i]<pred[j] else .5 if pred[i]==pred[j] else 0. for i,j in reliable)
            records.append(dict(arm=arm,receiver=rec,correct=wins,pairs=len(reliable),
                rho=spearmanr(pred[exact],true).statistic,tau_b=kendalltau(pred[exact],true).statistic,
                weak10_minus_strong8=pred[6]-pred[4],weak11_minus_strong9=pred[7]-pred[5]))
    reference={(r['arm'],r['receiver']):r for r in rows(submitted/'assay.csv')}
    for r in records:
        original=reference.pop((r['arm'],r['receiver']))
        for k,other in (('correct','correct'),('pairs','denominator'),('rho','rho'),('tau_b','tau_b'),
                        ('weak10_minus_strong8','cliff8_10'),('weak11_minus_strong9','cliff9_11')):
            np.testing.assert_allclose(r[k],float(original[other]),rtol=0,atol=1e-10)
    assert not reference
    table(out/'assay.csv',records)
    return dict(pools=112,nominal_pairs=28,threefold_pairs=len(reliable),point_values=len(exact),summaries=len(records))


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('group',choices=('pose','assay'))
    for key in ('work','out'):
        a.add_argument(key,type=Path)
    a.add_argument('--pose-results',type=Path)
    a.add_argument('--frames',type=Path)
    a.add_argument('--assay-results',type=Path,default=Path(__file__).resolve().parents[1]/'fkbp/results')
    args=a.parse_args()
    t,c=time.perf_counter(),time.process_time()
    r=(pose(args.work,args.frames,args.pose_results,args.out) if args.group=='pose' else assay(args.work,args.out,args.assay_results))
    r.update(wall_s=time.perf_counter()-t,cpu_s=time.process_time()-c)
    (args.out/(args.group+'_recount.json')).write_text(json.dumps(r,indent=2)+'\n')
