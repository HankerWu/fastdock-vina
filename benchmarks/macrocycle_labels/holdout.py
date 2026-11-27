"""Independent, graph-derived ligand start for human HSP90-alpha, PDB 1YET."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
import run_labels as trial

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path)
    for name in ['official','parent','fixed','full']:
        p.add_argument('--'+name, type=Path)
    p.add_argument('--prepare', action='store_true')
    a=p.parse_args()
    if a.prepare:
        sys.path.insert(0,str(trial.ROOT/'benchmarks/pair_gradient/pose'))
        import prepare
        prepare.CASES['1YET']=('HSP90-alpha','GDM','A')
        os.environ.setdefault('OPENMM_CPU_THREADS','1')
        prepare.prepare(a.work,'1YET')
        return
    src=a.work/'1YET'
    meta=json.loads((src/'preparation.json').read_text())
    if meta['receptor_status'] != 'ready':
        raise ValueError(meta['receptor_status'])
    raw=a.work/'docked'
    raw.mkdir(exist_ok=True)
    box=raw/'box.txt'
    box.write_text(''.join(f'center_{dim} = {x}\nsize_{dim} = 24\n'
                            for dim,x in zip('xyz',meta['center'])))
    results=[]
    for label in (0,3):
        lig=raw/f'GDM-{label}.pdbqt'
        trial.relabel(src/'ligand.pdbqt',lig,label)
        if f'CG{label}' not in lig.read_text():
            raise ValueError('Missing ring closure label')
        for engine in ('parent','fixed','full','official'):
            prefix=raw/f'{engine}-{label}'
            args=[str(getattr(a,engine).resolve()),'--config',str(box),'--ligand',str(lig),
                  '--receptor',str(src/'receptor.pdbqt'),'--cpu','1','--seed',str(trial.SEED),
                  '--max_evals',str(trial.CAP),'--num_modes','5','--energy_range','100',
                  '--exhaustiveness','1','--out',str(prefix.with_suffix('.pdbqt'))]
            if engine != 'official':
                args+=['--search_counters']
            start=time.perf_counter()
            log=subprocess.check_output(args,text=True,stderr=subprocess.STDOUT)
            wall=time.perf_counter()-start
            prefix.with_suffix('.log').write_text(log)
            poses=prefix.with_suffix('.pdbqt').read_text()
            xyz='\n'.join(line[:54] for line in poses.splitlines() if line.startswith(('ATOM','HETATM')))
            from geometry import coord_blocks
            modes = []
            for pose in coord_blocks(prefix.with_suffix('.pdbqt')):
                pos=np.array([pose[i][1] for i in meta['serial']])
                rmsd=np.sqrt(np.mean(np.sum((pos-np.array(meta['mapped']))**2,axis=-1),axis=1)).min()
                modes.append(round(float(rmsd),6))
            counts=trial.TRAJECTORY.search(log)
            refine=trial.REFINE.search(log)
            results.append(dict(engine=engine,label=label,top1=modes[0],best5=min(modes[:5]),modes=modes,
                                coordinate_hash=hashlib.sha256(xyz.encode()).hexdigest(),wall_s=round(wall,4),
                                hunt_evals=counts[2] if counts else None,second_evals=counts[4] if counts else None,
                                explicit_evals=refine[2] if refine else None,hunt_calls=counts[1] if counts else None,
                                second_calls=counts[3] if counts else None,explicit_calls=refine[1] if refine else None))
    for engine in ('fixed','full'):
        hashes={row['coordinate_hash'] for row in results if row['engine']==engine}
        if len(hashes)!=1:
            raise ValueError(f'{engine} not label invariant')
    with (a.work/'hsp90.csv').open('w') as handle:
        writer=csv.DictWriter(handle,results[0])
        writer.writeheader()
        writer.writerows(results)
if __name__=='__main__':
    main()
