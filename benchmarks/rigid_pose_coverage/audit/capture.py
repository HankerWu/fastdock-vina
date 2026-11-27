"""Follow native type-selected (blind) and earliest-near (labelled) forks from three frozen sources."""
import argparse
import json
import os
import resource
import subprocess
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

D = Path(__file__).resolve().parent
BENCH = D.parents[1]
sys.path[:0] = [str(BENCH/'fragment_anchor/audit'),str(BENCH/'cavity_triplet')]
from edge_audit import frames, equivalent_targets, raw_source
from forks import read_pool
from pose_runs import invoke, table, digest, counts
from pose_metrics import inspect_fork

CASES = ('3FMZ-A_93241','3FMZ-HA_93241','5TBM_93241')


def blind(cfg, work, case, seed, observer, costs):
    stem = case+'_coverage'
    out = work/'raw'/(stem+'.pdbqt')
    dump = work/'dump'/stem
    cmd = [sys.executable,D.parent/'engine/live.py','--engine',observer,'--arm','coverage','--dump',dump,'--',
           '--receptor',cfg['receptor'],'--ligand',cfg['ligand'],'--ligand_bonds',cfg['bonds'],'--cpu',1,
           '--seed',seed,'--exhaustiveness',2,'--max_evals',25000,'--num_modes',9,'--energy_range',100,
           '--pose_selection','late','--search_counters','--out',out]
    for axis,val in zip('xyz',cfg['center'],strict=True):
        cmd.extend(['--center_'+axis,val,'--size_'+axis,cfg['size']])
    log = out.with_suffix('.log')
    env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    wall, before = time.perf_counter(),resource.getrusage(resource.RUSAGE_CHILDREN)
    with log.open('w') as stream:
        subprocess.run(list(map(str,cmd)),stdout=stream,stderr=subprocess.STDOUT,env=env,check=True)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    costs.append(dict(stem=stem,wall_s=time.perf_counter()-wall,
        cpu_s=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,pool_sha256=digest(out),**counts(log)))
    for ext, path in (('.stages','stages.xyz'),('.inject','inject.log'),('.history','history')):
        out.with_suffix(ext).symlink_to(Path('..')/'dump'/stem/path)
    return stem


def main(parent, work, observer, result):
    raw = work/'raw'
    raw.mkdir(parents=True,exist_ok=False)
    (work/'dump').mkdir()
    cfgs = json.loads((parent/'inputs.json').read_text())['cells']
    costs, rows = [], []
    for case in CASES:
        name, seed = case.rsplit('_',1)
        cfg = cfgs[name]
        for key in ('receptor','bonds','ligand'):
            cfg[key] = parent/cfg[key]
        blind(cfg,work,case,seed,observer,costs)
    data = frames(work/'references')
    for case in CASES:
        name, seed = case.rsplit('_',1)
        cfg = cfgs[name]
        order, refs = equivalent_targets(name,cfg['ligand'],data)
        stem = case+'_coverage'
        bank_dir = work/'dump'/stem
        for source_file in ('source.xyz','source.q'):
            assert (bank_dir/source_file).read_bytes() == (parent/'dump'/(case+'_whole')/source_file).read_bytes()
        ids, x = raw_source(bank_dir/'source.xyz')
        targets = refs[:,[order.index(i) for i in ids]]
        with np.load(bank_dir/'bank.npz') as b:
            rot,center,triangle = b['matched_rotation'],b['center'],b['triangle_id']
        coords = (x-x.mean(0))@rot.transpose(0,2,1)+center[:,None,:]
        err = np.sqrt(np.minimum.reduce([np.square(coords-y).sum(2).mean(1) for y in targets]))
        scores = np.loadtxt(bank_dir/'scores.csv',delimiter=',',skiprows=1)[:,1]
        rank = np.argsort(scores,kind='stable')
        hits = np.flatnonzero(err[rank]<=2.)
        forks = [(stem,0)]
        if len(hits):
            j = rank[hits[0]]
            q = np.loadtxt(bank_dir/'source.q')[2:]
            orient = Rotation.from_matrix(rot[j])
            trans = orient.apply(q[:3]-x.mean(0))+center[j]
            quat = (orient*Rotation.from_quat(q[[4,5,6,3]])).as_quat()[[3,0,1,2]]
            labelled = stem+'_labelled'
            base = raw/labelled
            spec = base.with_suffix('.input')
            spec.write_text('1 1 '+' '.join(format(v,'.17g') for v in np.r_[trans,quat])+'\n')
            env=dict(FD_RIGID_INTERVENTION=spec, FD_RIGID_LOG=base.with_suffix('.inject'),
                     FD_RIGID_COORDS=base.with_suffix('.stages'),FD_RIGID_HISTORY=base.with_suffix('.history'))
            invoke(work,costs,labelled,observer,cfg,seed,env=env)
            forks.append((labelled,int(hits[0])))
        for tag, pos in forks:
            j = rank[pos]
            stages, history, lineage = inspect_fork(work,tag,dict(refs=refs,serials=order),1,1)
            np.testing.assert_allclose(err[j],lineage['selected_A'],rtol=0,atol=1e-10)
            pool = read_pool(raw/(tag+'.pdbqt'),order,refs)
            rows.append(dict(case=case,arm='labelled' if tag!=stem else 'blind',near=int(len(hits)),
                bank_min_A=err.min(),first_near_field_rank=int(hits[0]+1) if len(hits) else '',
                field_rank=pos+1,bank_index=j,triangle_id=triangle[j],score=scores[j],
                second_A=next((s['rmsd_A'] for s in stages if s['stage']=='second'),''),
                later_first_near_offer=next((r['step'] for r in history if r['step']>1 and r['offered_A']<=2.),''),
                top1_A=pool[0],oracle_A=min(pool),regret_A=pool[0]-min(pool),**lineage))
    result.mkdir(exist_ok=True,parents=True)
    table(result/'capture.csv',rows)
    table(result/'capture_cost.csv',costs)
    with tarfile.open(result/'capture_raw.tar.xz','w:xz') as archive:
        for sub in ('dump','raw'):
            archive.add(work/sub,arcname=sub)
    print(digest(result/'capture_raw.tar.xz'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('parent','work','observer','result'):
        parser.add_argument(name,type=lambda s:Path(s).resolve())
    args=parser.parse_args()
    main(args.parent,args.work,args.observer,args.result)
