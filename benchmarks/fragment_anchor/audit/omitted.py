"""After prediction, quantify how different the discarded near poses are from retained candidates."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import hypergeom

from edge_audit import Pocket, align, raw_source, table
from scipy.spatial.distance import pdist


def run(work, out):
    cells=json.loads((work/'inputs.json').read_text())['cells']
    cases=list(csv.DictReader((out/'gates_pose.csv').open()))
    details=[]
    with np.load(out/'pose_edges.npz') as curves:
        for row in cases:
            case,arm=row['case'],row['arm']
            key=case+'_'+arm
            e=curves[key]
            near=np.flatnonzero(e<=2.)
            assert len(near)==int(row['before_le2']) and sum(near<1024)==int(row['after_le2'])
            lost=near[near>=1024]
            entry=dict(case=case,arm=arm,all_le2=len(near),retained_le2=int(sum(near<1024)),omitted=len(lost),
                       uniform_1024_hit_probability=hypergeom.sf(0,len(e),len(near),1024),
                       omitted_to_retained_min_A='',omitted_to_retained_median_A='',closest_retained_crystal_A='')
            if len(lost):
                c=cells[case.split('_')[0]]
                pocket=Pocket.from_pdbqt(work/c['receptor'],c['center'],c['size'])
                path=work/'dump'/key
                ids,x=raw_source(path/'source.xyz')
                with np.load(path/'bank.npz') as archive:
                    tri=archive['ligand_triangle']
                    kept_centers,kept_rot=archive['center'],archive['matched_rotation'].transpose(0,2,1)
                a=x[tri]
                rank=np.argsort(np.round(((pocket.edges-pdist(a))**2).sum(1),10),kind='stable')
                y=pocket.sites[pocket.triangles[rank[lost]]]
                allrot=align(a,y)
                allcenters=(x.mean(0)-a.mean(0))@allrot+y.mean(1)
                centered=x-x.mean(0)
                gram=centered.T@centered/len(x)
                nearest,neighbours=[],[]
                for r,t in zip(allrot,allcenters,strict=True):
                    d=kept_rot-r
                    sq=(np.square(kept_centers-t).sum(1)+np.einsum('nai,ab,nbi->n',d,gram,d))
                    nn=sq.argmin()
                    nearest.append(np.sqrt(sq[nn]))
                    neighbours.append(e[nn])
                entry.update(omitted_to_retained_min_A=min(nearest),omitted_to_retained_median_A=np.median(nearest),
                             closest_retained_crystal_A=neighbours[np.argmin(nearest)])
            details.append(entry)
    table(out/'omitted.csv',details)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path)
    p.add_argument('out',type=Path)
    a=p.parse_args()
    run(a.work,a.out)
