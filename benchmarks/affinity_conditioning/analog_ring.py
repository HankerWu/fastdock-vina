"""Compare the shared pipecolate ring of prepared 5/6 to the 8/9 analog crystallographic rings."""
import argparse
import csv
import hashlib
import itertools
import json
import time
from pathlib import Path

import numpy as np
from rdkit import Chem

from verify import BENCH, Ligand, members


def loop(cif, name):
    block, = [b for b in cif.split("loop_") if b.lstrip().startswith(f"_chem_comp_{name}.")]
    lines = block[:block.index("#")].strip().splitlines()
    names = [s.strip().split(".")[1] for s in lines if s.startswith("_")]
    return [dict(zip(names,s.split(),strict=True)) for s in lines if not s.startswith("_")]


def crystal(files, pdb, comp):
    text = files[f"raw/{pdb}.pdb"]
    cif = files[f"raw/{comp}.cif"]
    atoms = {r["atom_id"]:r for r in loop(cif,"atom")}
    bonds = {frozenset((r["atom_id_1"],r["atom_id_2"])):r["value_order"] for r in loop(cif,"bond")}
    ring = ["N7","C6","C5","C4","C3","C2"]
    edges = [*zip(ring, ring[1:]+ring[:1]), ("C2","C1"),("N7","C8"),("C1","O2")]
    if any(bonds[frozenset(e)]!="SING" for e in edges):
        raise ValueError("Unexpected pipecolate connectivity")
    if atoms["C2"]["pdbx_stereo_config"]!="S" or bonds[frozenset(("C1","O1"))]!="DOUB":
        raise ValueError("Invalid alpha-carbon or ester")
    if bonds[frozenset(("C8","O3"))]!="DOUB":
        raise ValueError("Missing N-acyl carbonyl")
    xyz = {}
    for row in text.splitlines():
        if row.startswith("HETATM") and row[17:20]==comp and row[21]=="A" and int(row[22:26])==108:
            name=row[12:16].strip()
            if name in xyz or row[16]!=" " or float(row[54:60])!=1.0:
                raise ValueError("Unexpected ambiguous ligand coordinate")
            xyz[name]=np.array([float(row[k:k+8]) for k in (30,38,46)])
    return ring, np.array([xyz[name] for name in ring]), xyz["C1"]


def fitted_rmsds(source, reference):
    x=source-source.mean(axis=0)
    y=reference-reference.mean(axis=0)
    h=x.T@y
    u, _, vt=np.linalg.svd(h)
    u[:,-1]*=np.linalg.det(u@vt)
    rmsd=np.linalg.norm(x@(u@vt)-y)/np.sqrt(len(x))
    # Maximize the alignment trace independently via Horn's unit quaternion.
    tr=np.trace(h)
    antisymmetric=np.array([h[1,2]-h[2,1],h[2,0]-h[0,2],h[0,1]-h[1,0]])
    davenport=np.block([[np.array([[tr]]),antisymmetric[None,:]],
                        [antisymmetric[:,None],h+h.T-tr*np.eye(3)]])
    largest=np.linalg.eigvalsh(davenport)[-1]
    quaternion_rmsd=np.sqrt(max(0., (np.sum(x*x)+np.sum(y*y)-2*largest)/len(x)))
    return rmsd,quaternion_rmsd


def run(out):
    clock=time.perf_counter(),time.process_time()
    source=BENCH/"cubic_grid/human_panel"
    inputs=members(source/"depositions.tar.xz",["raw/1FKG.pdb","raw/1FKH.pdb","raw/SB3.cif","raw/SBX.cif"])
    ligands=members(source/"affinity-pools.tar.xz",["affinity/5.pdbqt","affinity/6.pdbqt"])
    structures={}
    pattern=Chem.MolFromSmarts("[N;R]1[C;R][C;R][C;R][C;R][C;R]1C(=O)O")
    for n in (5,6):
        ligand=Ligand(ligands[f"affinity/{n}.pdbqt"])
        match,=ligand.graph.GetSubstructMatches(pattern)
        ids=list(match[:6])
        structures[f"prepared_{n}"]=([ligand.order[i] for i in ids],
                                     ligand.reference[ids],ligand.reference[match[6]])
    for n,pdb,comp in ((8,"1FKG","SB3"),(9,"1FKH","SBX")):
        structures[f"crystal_{n}_{pdb}"]=crystal(inputs,pdb,comp)
    signs={}
    for name,(_,xyz,outer) in structures.items():
        # Ordered alpha-carbon heavy neighbors are ester carbonyl, ring C3, ring N.
        signs[name]=float(np.linalg.det(np.stack((outer,xyz[4],xyz[0]))-xyz[5]))
    if len(set(np.sign(v) for v in signs.values()))!=1 or min(map(abs,signs.values()))<0.5:
        raise ValueError("Deposited and prepared alpha-center configurations disagree")
    out.mkdir(parents=True,exist_ok=True)
    results=[]
    for first,second in itertools.combinations(structures,2):
        r1,r2=fitted_rmsds(structures[first][1], structures[second][1])
        if abs(r1-r2)>1e-10:
            raise ValueError("Independent proper-rigid fits disagree")
        results.append(dict(first=first, second=second, svd_rmsd_A=r1, quaternion_rmsd_A=r2))
    with (out/"analog_ring.csv").open("w", newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=list(results[0]),lineterminator="\n")
        w.writeheader()
        w.writerows(results)
    records=dict(source_sha256={k:hashlib.sha256(v.encode()).hexdigest() for k,v in {**inputs,**ligands}.items()},
                 rings={k:dict(atom_sequence=ids, xyz_A=xyz.tolist(), alpha_volume_A3=signs[k])
                        for k,(ids,xyz,_) in structures.items()},
                 wall_s=time.perf_counter()-clock[0],cpu_s=time.process_time()-clock[1])
    (out/"analog_ring_source.json").write_text(json.dumps(records,indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    run(p.parse_args().out)
