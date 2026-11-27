"""Check dependence of each local diagnostic on enlargement of the original docking box."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import numpy as np
from swap import coords, energy, rmsd, write


def check(work,binary):
    folder=work/"box24"
    folder.mkdir()
    with (work/"factorial.csv").open() as stream:
        rows=list(csv.DictReader(stream))
    stats={}
    output=[]
    for row in rows:
        pdb=row["pdb"]
        key=row["artifact"]
        if key not in stats:
            site=work/"source"/pdb
            meta=json.loads((site/"preparation.json").read_text())
            pose=folder/(key+".pdbqt")
            cmd=[str(binary),"--cpu","1","--seed","43103","--local_only","--receptor",str(site/"receptor.pdbqt"),
                 "--ligand",str(work/"runs"/(key+".pdbqt")),"--out",str(pose)]
            for axis,v in zip("xyz",meta["center"]):
                cmd.extend(["--center_"+axis,str(v),"--size_"+axis,"24"])
            proc=subprocess.run(cmd,capture_output=True,text=True)
            (folder/(key+".log")).write_text(proc.stdout+proc.stderr)
            if proc.returncode:
                assert "outside the grid box" in proc.stderr,proc.stderr
                stats[key]=dict(status="outside",final_rmsd24="",local_energy24="")
            else:
                points=coords(pose.read_text())
                x=np.array([points[i] for i in meta["serial"]])
                stats[key]=dict(status="ok",final_rmsd24=rmsd(x,np.array(meta["mapped"])),local_energy24=energy(proc.stdout))
        output.append({**{k:row[k] for k in ("pdb","cap","seed","rank","selection","arm","artifact","final_rmsd")},**stats[key]})
    write(work/"box24.csv",output)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("binary",type=Path)
    a=p.parse_args()
    check(a.work.resolve(),a.binary.resolve())
