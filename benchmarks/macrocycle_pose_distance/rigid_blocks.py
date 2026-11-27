"""Check whether exact deposited positions are attainable by this prepared torsion tree."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"macrocycle_proposal"/"external"))
from external_score import poses, write


def blocks(path):
    tree = []
    rigid = defaultdict(list)
    for line in path.read_text().splitlines():
        if line.startswith("BRANCH "):
            tree.append(int(line.split()[2]))
        elif line.startswith("ENDBRANCH "):
            tree.pop()
        elif line.startswith(("ATOM  ","HETATM")):
            rigid[tuple(tree)].append(int(line[6:11]))
    return rigid


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    rows=[]
    for pdb in ("2J9M","3VHA","4EGK"):
        meta=json.loads((a.inputs/pdb/"preparation.json").read_text())
        path=a.inputs/pdb/"ligand.pdbqt"
        groups=blocks(path)
        coords=next(poses(path))
        serial=meta["serial"]
        for k,ref in enumerate(np.asarray(meta["mapped"])):
            pairs=[]
            for group in groups.values():
                chemical=[i for i in group if i in serial]
                for h,i in enumerate(chemical):
                    for j in chemical[:h]:
                        moving=np.linalg.norm(coords[i][1]-coords[j][1])
                        fixed=np.linalg.norm(ref[serial.index(i)]-ref[serial.index(j)])
                        pairs.append((abs(moving-fixed),i,j,moving,fixed))
            gap,i,j,prepared,deposited=max(pairs)
            rows.append(dict(pdb=pdb,mapping=k+1,serial_i=i,serial_j=j,prepared_A=prepared,
                             deposited_A=deposited,max_invariant_distance_difference_A=gap))
    write(a.out,rows)
