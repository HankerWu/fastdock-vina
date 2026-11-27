"""Measure geometry within query PDBQT rigid fragments and MALT1's fused aromatic butterfly."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.spatial.distance import pdist,squareform
from rdkit import Chem
from swap import coords, transform, write


def fragments(template,serials):
    groups={}
    stack=[]
    for line in template.splitlines():
        if line=="ROOT":
            stack=[0]
        elif line.startswith("BRANCH"):
            stack.append(int(line.split()[2]))
        elif line.startswith("ENDBRANCH"):
            stack.pop()
        elif line.startswith(("ATOM  ","HETATM")):
            s=int(line[6:11])
            if s in serials:
                groups.setdefault(stack[-1],[]).append(serials.index(s))
    return groups


def main(work):
    summary=[]
    rings=[]
    for pdb in ("2YIY","4I1R"):
        site=work/"source"/pdb
        meta=json.loads((site/"preparation.json").read_text())
        tpl=(site/"ligand.pdbqt").read_text()
        gr=fragments(tpl,meta["serial"])
        raw=coords((work/"source/crystal_diagnostic"/f"{pdb}-{pdb}-best_projected-input.pdbqt").read_text())
        x=np.array([raw[s] for s in meta["serial"]])
        y=min(np.array(meta["mapped"]),key=lambda m:((x-m)**2).sum())
        for axis,ids in gr.items():
            if len(ids)<3:
                continue
            a,b=x[ids],y[ids]
            rot,t=transform(a,b)
            delta=abs(pdist(a)-pdist(b))
            i,j=np.unravel_index(squareform(delta).argmax(),(len(ids),len(ids)))
            summary.append(dict(pdb=pdb,group=axis,atoms=len(ids),fitted_rmsd_A=float(np.sqrt(((a@rot+t-b)**2).sum(1).mean())),
                                max_distance_diff_A=max(delta),pdbqt_serial_a=meta["serial"][ids[i]],pdbqt_serial_b=meta["serial"][ids[j]],
                                projected_length_A=np.linalg.norm(a[i]-a[j]),crystal_length_A=np.linalg.norm(b[i]-b[j])))
        if pdb!="4I1R":
            continue
        m=Chem.MolFromSmiles(meta["canonical_smiles"])
        aromatic=[list(r) for r in m.GetRingInfo().AtomRings() if len(r)==6 and all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
        assert len(aromatic)==2 and any(set(aromatic[0]+aromatic[1])<=set(ids) for ids in gr.values())
        # Use absolute normals: the acute dihedral describes bending magnitude, not which way the ring bends.
        for name,positions in (("projected",x),("crystal",y)):
            normals=[np.linalg.svd(positions[ids]-positions[ids].mean(axis=0))[2][-1] for ids in aromatic]
            angle=np.degrees(np.arccos(np.clip(abs(normals[0]@normals[1]),0,1)))
            rings.append(dict(pdb=pdb,geometry=name,plane_angle_deg=angle,ring_serials=str([[meta["serial"][i] for i in ring] for ring in aromatic])))
    write(work/"rigid_fragments.csv",summary)
    write(work/"malt1_aromatics.csv",rings)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    main(p.parse_args().work)
