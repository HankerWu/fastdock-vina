"""Score every mode of the original G0 output on identical coordinates and topology."""
import argparse
import csv
from pathlib import Path
import re
import subprocess
import tempfile

from run_labels import relabel

def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('parent','fixed','full','bace-work','hsp-work','upstream-data'):
        p.add_argument('--'+name, required=True, type=Path)
    a=p.parse_args()
    data=a.upstream_data/'solution'
    cases=[
        ('5QCT',a.bace_work/'raw/single-parent-0.pdbqt',data/'BACE_1_ligand.pdbqt',
         data/'BACE_1_receptor.pdbqt',data/'BACE_1_receptor_vina_box.txt'),
        ('1YET',a.hsp_work/'docked/parent-0.pdbqt',a.hsp_work/'1YET/ligand.pdbqt',
         a.hsp_work/'1YET/receptor.pdbqt',a.hsp_work/'docked/box.txt')]
    rows=[]
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp)
        for pdb,pool,lig,rec,box in cases:
            for rank,pose in enumerate(re.findall(r'^MODEL .*?^ENDMDL',pool.read_text(),re.M|re.S),1):
                original=tmp/'original.pdbqt'
                original.write_text(pose+'\n')
                for label in (0,3):
                    prepared=tmp/'prepared.pdbqt'
                    coordinates=tmp/'coordinates.pdbqt'
                    relabel(lig,prepared,label)
                    relabel(original,coordinates,label)
                    for engine in ('parent','fixed','full'):
                        args=[str(getattr(a,engine).resolve()),'--score_only','--ligand',str(prepared),
                              '--pose_coordinates',str(coordinates),'--unbound_energy','0','--cpu','1',
                              '--seed','458921','--config',str(box),'--receptor',str(rec)]
                        log=subprocess.check_output(args,text=True,stderr=subprocess.STDOUT)
                        value=re.search(r'Estimated Free Energy of Binding\s*:\s*([-\d.]+)',log)[1]
                        rows.append(dict(pdb=pdb,rank=rank,label=label,engine=engine,energy=float(value)))
    with (a.bace_work/'fixed_pool.csv').open('w') as handle:
        writer=csv.DictWriter(handle,rows[0])
        writer.writeheader()
        writer.writerows(rows)

if __name__ == '__main__':
    main()
