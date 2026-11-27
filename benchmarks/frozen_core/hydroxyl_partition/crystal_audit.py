"""Independent ESR2 reference construction from the deposited PDB ligand records.

Infer the observed ligand's neighbor graph directly from coordinates, assign
bond orders from SMILES, and enumerate stereo-preserving exact graph mappings.
The production preparation's CCD-name-to-coordinate map is never read.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import tarfile
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem

from scorer import atoms


def assess():
    here = Path(__file__).resolve().parent
    with (here/"results/modes.csv").open() as f:
        old = {(r["query"],r["cap"],r["seed"],int(r["rank"])):float(r["rmsd_A"])
               for r in csv.DictReader(f) if r["panel"]=="esr2"}
    reports=[]
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        with tarfile.open(here/"esr2.tar.xz") as t:
            t.extractall(root,filter="data")
        for code,component in (("1QKM","GEN"),("5TOA","EST")):
            text=(root/code/"ligand.pdbqt").read_text()
            smi = next(l.split()[2] for l in text.splitlines() if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX "))
            query=Chem.MolFromSmiles(smi)
            serial={}
            for l in text.splitlines():
                if l.startswith("REMARK SMILES IDX "):
                    indices=list(map(int,l.split()[3:]))
                    serial.update(zip(indices[::2],indices[1::2],strict=True))
            lines=[l for l in (root/"raw"/(code+".pdb")).read_text().splitlines() if l.startswith("HETATM")
                   and l[21]=="A" and l[17:20]==component and l[16] in (" ","A") and l[76:78].strip() not in ("H","D")]
            observed=Chem.MolFromPDBBlock("\n".join(lines)+"\nEND\n",sanitize=False,removeHs=False)
            assert observed is not None and observed.GetNumAtoms()==query.GetNumAtoms()
            observed=AllChem.AssignBondOrdersFromTemplate(query,observed)
            Chem.AssignStereochemistryFrom3D(observed,replaceExistingTags=True)
            maps=observed.GetSubstructMatches(query,useChirality=True,uniquify=False,maxMatches=1001)
            assert maps and len(maps)<1001
            reference=observed.GetConformer().GetPositions()
            delta=[]
            for pool in sorted((root/code).glob("*-*-*.pdbqt")):
                _,cap,seed=pool.stem.split("-")
                for i,block in enumerate((b for b in pool.read_text().split("ENDMDL") if b.strip()),1):
                    positions=atoms(block)
                    coords=[positions[serial[k+1]].xyz for k in range(query.GetNumAtoms())]
                    rms=min(math.sqrt(math.fsum(math.dist(c,reference[j])**2 for c,j in zip(coords,m,strict=True))/len(m)) for m in maps)
                    delta.append(abs(rms-old[code,cap,seed,i]))
            reports.append(dict(pdb=code,heavy=query.GetNumAtoms(),stereo_mappings=len(maps),poses=len(delta),max_rmsd_diff_A=max(delta)))
    assert sum(r["poses"] for r in reports)==len(old)
    assert all(r["max_rmsd_diff_A"]<1e-12 for r in reports)
    return reports


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    args=p.parse_args()
    args.out.write_text(json.dumps(assess(),indent=2)+"\n")
