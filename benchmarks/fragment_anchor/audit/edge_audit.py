"""Evaluation-only geometry in the receiver frame, including rejected cavity correspondences."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile
import time

import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
TRIPLET = BENCH / 'cavity_triplet'
sys.path[:0] = [str(HERE.parent / 'engine'), str(BENCH / 'cavity_field_rank/audit')]
from anchors import select_anchor
from direct_space import PDB, CODES, deposited_frames, equivalent_targets, reference, table
from pocket_triangle import Pocket, _ligand_triangle


def frames(work):
    data = deposited_frames(work)
    fabp = {'2NNQ':'T4B', '3FR2':'8CA'}
    with tarfile.open(BENCH / 'cubic_grid/human_panel/depositions.tar.xz') as source:
        for pdb, component in fabp.items():
            for entry in (pdb+'.pdb',component+'.cif'):
                source.extract('raw/'+entry, work, filter='data')
            data[pdb] = reference(work,pdb,component)
    PDB.update({'2NNQ':('2NNQ','3FR2'),'3FR2':('3FR2','2NNQ')})
    return data


def raw_source(path):
    header, line = path.read_text().splitlines()
    assert line.split()[:3] == ['source','1','1']
    ids = np.array(header.split()[1:],int)
    return np.sort(ids), np.array(line.split()[3:],float).reshape(-1,3)[np.argsort(ids)]


def align(x, y):
    """Row-vector proper rotations fitted to landmark correspondences; y is batched."""
    source = x-x.mean(0)
    dest = y-y.mean(axis=1,keepdims=True)
    u,_,vt = np.linalg.svd(np.einsum('ji,njk->nik',source,dest))
    u[:,:,-1] *= np.sign(np.linalg.det(u@vt))[:,None]
    return u@vt


def fits(x, ref, anchor):
    """Rigidly place the actual source on ideal *labelled* crystal landmarks."""
    estimates = []
    for target in ref:
        a,b = x[anchor],target[anchor]
        rot,_ = Rotation.align_vectors(b-b.mean(0),a-a.mean(0))
        actual = rot.apply(x-a.mean(0)) + b.mean(0)
        estimates.append(min(np.linalg.norm(actual-y)/np.sqrt(len(x)) for y in ref))
    return min(estimates)


def measure(case, arm, x, ids, ref, pocket, lig, bonds, check=None):
    if arm=='whole':
        tri = _ligand_triangle(x)
        comp, fallback = ids,False
    else:
        tri, comp, fallback = select_anchor(x,ids,lig,bonds)
    a=x[tri]
    error2 = ((pocket.edges-pdist(a))**2).sum(1)
    ix = np.argsort(np.round(error2,10),kind='stable')
    triples = pocket.triangles[ix]
    # Fixed receiver frame, optimize only chemical atom permutations. Never align predicted pose.
    errors = np.empty(len(ix))
    for start in range(0,len(ix),128):
        sub = slice(start,start+128)
        y = pocket.sites[triples[sub]]
        rotation = align(a,y)
        posed = (x-a.mean(0))@rotation + y.mean(axis=1)[:,None,:]
        errors[sub] = np.sqrt(np.min([np.square(posed-t).sum(2).mean(1) for t in ref],axis=0))
    if check is not None:
        assert np.array_equal(triples[:1024],check['site_triangles'])
        assert np.array_equal(ids[tri],check['anchor_serials'])
        xyz = (x-x.mean(0))@check['matched_rotation'].transpose(0,2,1)+check['center'][:,None,:]
        err2 = np.sqrt(np.min([np.square(xyz-t).sum(2).mean(1) for t in ref],axis=0))
        np.testing.assert_allclose(errors[:1024],err2,rtol=0,atol=1e-9)
    before = np.flatnonzero(errors<=2)
    after = errors[:1024]
    row = dict(case=case,arm=arm,n_sites=len(pocket.sites),correspondences=len(ix),
               anchor_serials=' '.join(map(str,ids[tri])),component_size=len(comp),fallback=int(fallback),
               diameter_A=np.linalg.norm(a[0]-a[1]),double_area_A2=np.linalg.norm(np.cross(a[1]-a[0],a[2]-a[0])),
               ideal_rigid_A=fits(x,ref,np.arange(len(x))),ideal_anchor_on_crystal_A=fits(x,ref,tri),
               before_min_A=min(errors),before_le2=len(before),before_le3=int(sum(errors<=3)),
               best_edge_rank=int(errors.argmin()+1),first_le2_edge_rank=(int(before[0]+1) if len(before) else ''),
               after_min_A=min(after),after_le2=int(sum(after<=2)),after_le3=int(sum(after<=3)),
               after_best_id=int(after.argmin()))
    if row['before_min_A']+1e-8 < row['ideal_rigid_A']:
        raise ValueError('impossible lower bound')
    return row, errors


def pose_cases(p, work, out):
    data = frames(work)
    manifest = json.loads((p/'inputs.json').read_text())['cells']
    cases = sorted(d.parent.name.rsplit('_',1)[0] for d in (p/'dump').glob('*_whole/bank.npz'))
    assert len(cases)==16
    records,raw=[],{}
    for case in cases:
        name,seed=case.split('_')
        c=manifest[name]
        for file,sha in c['sha256'].items():
            assert hashlib.sha256((p/file).read_bytes()).hexdigest()==sha
        ligand,bonds=p/c['ligand'],p/c['bonds']
        order,target=equivalent_targets(name,ligand,data)
        base=p/'dump'/f'{case}_whole'
        ids,x=raw_source(base/'source.xyz')
        ref=target[:,[order.index(s) for s in ids]]
        pocket=Pocket.from_pdbqt(p/c['receptor'],c['center'],c['size'])
        for arm in ('whole','fragment'):
            d=p/'dump'/f'{case}_{arm}'
            assert (base/'source.xyz').read_bytes()==(d/'source.xyz').read_bytes()
            with np.load(d/'bank.npz') as bank:
                report,err=measure(case,arm,x,ids,ref,pocket,ligand,bonds,bank)
            records.append(report)
            raw[f'{case}_{arm}']=err.astype('float32')
    old={r['case']:r for r in csv.DictReader((BENCH/'cavity_field_rank/audit/results/representation.csv').open())}
    for row in records:
        if row['arm']=='whole' and row['case'] in old:
            expect=old[row['case']]
            np.testing.assert_allclose([row['after_min_A'],row['ideal_rigid_A']],
                                       [float(expect['matched_min_A']),float(expect['ideal_rigid_A'])],atol=1e-8)
    table(out/'gates_pose.csv',records)
    np.savez_compressed(out/'pose_edges.npz',**raw)


def aurora(work,out):
    aur=BENCH/'receptor_linkage/aurora'
    with tarfile.open(aur/'references.tar.xz') as tar:
        tar.extractall(work,filter='data')
    for pkg, paths in [(aur/'docking.tar.xz', ['prepared/at9283.sdf','prepared/at9283.pdbqt']),
                       (aur/'probe_raw.tar.xz',['raw/at9283_chemical.bonds']),
                       (aur/'internal_loops/evidence.tar.xz',['template1OL7/trans_model/released/receptor.pdbqt'])]:
        with tarfile.open(pkg) as tar:
            for path in paths:
                tar.extract(path,work,filter='data')
    sys.path.insert(0,str(aur/'internal_loops'))
    import replay
    ligand,bonds,order,target=replay.parameters(work,'at9283')
    receptor=work/replay.PROTEINS['1OL7_released']
    c=json.loads((TRIPLET/'results/pose_inputs.json').read_text())['systems']['at9283']
    for k,path in [('ligand',ligand),('bonds',bonds),('receptor',receptor)]:
        assert hashlib.sha256(path.read_bytes()).hexdigest()==c[k]
    old={r['source']:r for r in csv.DictReader((TRIPLET/'results/posthoc_edge_gate.csv').open())}
    records,raw=[],{}
    pocket=Pocket.from_pdbqt(receptor,c['center'],c['size'])
    with tarfile.open(TRIPLET/'pose_raw.tar.xz') as tar:
        for member in tar:
            if not member.name.startswith('sources/at9283_61231_'):
                continue
            case=Path(member.name).stem
            with np.load(io.BytesIO(tar.extractfile(member).read())) as src:
                ix=np.argsort(src['serials'])
                ids,x=src['serials'][ix],src['xyz'][ix]
            ref=target[:,[list(order).index(s) for s in ids]]
            for arm in ('whole','fragment'):
                row,err=measure(case,arm,x,ids,ref,pocket,ligand,bonds)
                if arm=='whole':
                    historical=old[case]
                    for key,h in [('before_min_A','before_gate_min_A'),('after_min_A','after_gate_min_A')]:
                        np.testing.assert_allclose(row[key],float(historical[h]),atol=1e-9)
                records.append(row)
                raw[f'{case}_{arm}']=err.astype('float32')
    assert len(records)==36
    table(out/'gates_at9283.csv',records)
    np.savez_compressed(out/'at9283_edges.npz',**raw)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('group',choices=['pose','at9283'])
    parser.add_argument('work',type=Path)
    parser.add_argument('out',type=Path)
    parser.add_argument('--pose',type=Path)
    a=parser.parse_args()
    a.work.mkdir(parents=True,exist_ok=True)
    a.out.mkdir(parents=True,exist_ok=True)
    begin,cpu=time.perf_counter(),time.process_time()
    if a.group=='pose':
        pose_cases(a.pose,a.work,a.out)
    else:
        aurora(a.work,a.out)
    (a.out/(a.group+'_geometry_time.json')).write_text(json.dumps(
        dict(wall_s=time.perf_counter()-begin,cpu_s=time.process_time()-cpu),indent=2)+'\n')
