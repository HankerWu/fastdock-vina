"""Restore three deposited OLC molecules as rigid receptor components of human GLUT3."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import gemmi
import numpy as np
from meeko import MoleculePreparation
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[2]
sys.path.insert(0,str(BENCH/"hydration_transfer"))
from panel import cross, timed_process, write_csv

OLC_URL = "https://files.rcsb.org/ligands/download/OLC.cif"


def rigid_olc(ccd, records, seq, clock, dry_lines):
    block = gemmi.cif.read_file(str(ccd)).sole_block()
    atoms = block.get_mmcif_category("_chem_comp_atom.")
    bonds = block.get_mmcif_category("_chem_comp_bond.")
    g = Chem.RWMol()
    ids = {}
    for name,element,charge in zip(atoms["atom_id"],atoms["type_symbol"],atoms["charge"]):
        if element!="H":
            atom = Chem.Atom(element)
            atom.SetFormalCharge(int(charge))
            ids[name] = g.AddAtom(atom)
    for a,b,order in zip(bonds["atom_id_1"],bonds["atom_id_2"],bonds["value_order"]):
        if a in ids and b in ids:
            g.AddBond(ids[a],ids[b],{"SING":Chem.BondType.SINGLE,"DOUB":Chem.BondType.DOUBLE}[order])
    graph = g.GetMol()
    Chem.SanitizeMol(graph)
    ref = {l[12:16].strip():l for l in records}
    if len(ref)!=len(records) or set(ids)!=set(ref) or len(ids)!=25:
        raise ValueError("Incomplete or ambiguous deposited OLC")
    conf = Chem.Conformer(len(ids))
    for name,i in ids.items():
        conf.SetAtomPosition(i,cross.coords(ref[name]))
    graph.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(graph)
    desc = block.get_mmcif_category("_pdbx_chem_comp_descriptor.")
    expected = next(s for s,t in zip(desc["descriptor"],desc["type"]) if t=="SMILES_CANONICAL")
    computed=Chem.MolToSmiles(graph)
    expected=Chem.MolToSmiles(Chem.MolFromSmiles(expected))
    # Some protein crystal lipids have mis-modelled stereo; record as a caveat.
    full = Chem.AddHs(graph,addCoords=True)
    prop = AllChem.MMFFGetMoleculeProperties(full,mmffVariant="MMFF94s")
    field = AllChem.MMFFGetMoleculeForceField(full,prop)
    for i in ids.values():
        field.AddFixedPoint(i)
    status = field.Minimize(maxIts=2000)
    if status:
        raise ValueError(f"OLC {seq}: hydrogen relaxation failed")
    prep, = MoleculePreparation().prepare(full)
    labels = {i:n for n,i in ids.items()}
    lines = []
    positions = []
    for a in prep.atoms:
        if a.is_ignore:
            continue
        name = labels[a.index] if a.index in labels else "H"+labels[a.graph[0]].removeprefix("O")
        x,y,z = a.coord
        if a.index in labels and max(abs(np.asarray(cross.coords(ref[name]))-a.coord))>1.e-10:
            raise ValueError("Lipid heavy atom moved")
        clock += 1
        lines.append(f"HETATM{clock:5d} {name:>4} OLC A{seq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}    {a.charge:6.3f} {a.atom_type:<2}")
        if a.atomic_num!=1:
            positions.append(a.coord)
    if len(lines)!=27:
        raise ValueError("Expected 25 heavy atoms and two polar hydrogens")
    protein = [cross.coords(l) for l in dry_lines if l.startswith("ATOM") and l.split()[-1] not in {"HD","H"}]
    dist = cdist(positions,protein)
    return lines,clock,dict(residue=seq,complete_heavy_atoms=len(positions),polar_H=len(lines)-len(positions),
             canonical_smiles=expected,observed_3d_smiles=computed,stereo_matches_CCD=computed==expected,
             deposited_occupancies=sorted({float(l[54:60]) for l in records}),
             alternate_locations=sorted({l[16] for l in records}),shortest_protein_lipid_A=float(dist.min()))


def prepare(w,o):
    o.mkdir(exist_ok=True,parents=True)
    start=time.perf_counter()
    ccd=o/"OLC.cif"
    if not ccd.exists():
        ccd.write_bytes(gzip.decompress((HERE/"OLC.cif.gz").read_bytes()))
    ccd_sha=cross.fetch(ccd,OLC_URL)
    dry=(w/"4ZW9/receptor.pdbqt").read_text().splitlines()
    lines=(w/"raw/4ZW9.pdb").read_text().splitlines()
    clock=max(int(l[6:11]) for l in dry if l.startswith(("ATOM  ","HETATM")))
    all_lipids=[]
    receipts=[]
    for seq in (501,502,503):
        atoms=[l for l in lines if l.startswith("HETATM") and l[17:20]=="OLC" and l[21]=="A" and int(l[22:26])==seq]
        lipid,clock,summary=rigid_olc(ccd,atoms,seq,clock,dry)
        for code in ("GLC","BGC"):
            sugar=[l for l in lines if l.startswith("HETATM") and l[17:20]==code and l[21]=="A"]
            d=cdist([cross.coords(l) for l in atoms],[cross.coords(l) for l in sugar])
            i,j=np.unravel_index(np.argmin(d),d.shape)
            summary[code+"_nearest"]=[atoms[i][12:16].strip(),sugar[j][12:16].strip(),float(d[i,j])]
        all_lipids+=lipid
        receipts.append(summary)
    (o/"lipids.pdbqt").write_text("\n".join(all_lipids)+"\n")
    (o/"4ZW9-OLC.pdbqt").write_text("\n".join(dry+all_lipids)+"\n")
    write_csv(o/"preparation.csv",receipts)
    (o/"source.json").write_text(json.dumps(dict(source=OLC_URL,CCD_sha256=ccd_sha,seconds=time.perf_counter()-start),indent=2)+"\n")


def dock(w,o,binary):
    meta=json.loads((w/"metadata.json").read_text())
    (o/"runs").mkdir(exist_ok=True)
    c=[]
    for code,cap,seed in itertools.product(("GLC","BGC"),(3000,12000),(43103,43117)):
        m=meta["4ZW9-"+code]
        stem=f"4ZW9-{code}_{cap}_{seed}_olc"
        cmd=[str(binary),"--ligand",str(w/"4ZW9"/(code+".pdbqt")),"--receptor","4ZW9-OLC.pdbqt",
             "--seed",str(seed),"--cpu","1","--exhaustiveness","4","--max_evals",str(cap),
             "--num_modes","9","--energy_range","100","--out","runs/"+stem+".pdbqt","--search_counters"]
        for ch,x in zip("xyz",m["center"]):
            cmd.extend(["--center_"+ch,str(x),"--size_"+ch,"24"])
        log=o/"runs"/(stem+".log")
        timing=timed_process(cmd,o,log)
        trace=log.read_text()
        parts=re.findall(r'hunt_evals=(\d+).*?second_evals=(\d+)',trace)
        refined=re.findall(r'explicit_evals=(\d+)',trace)
        c.append(dict(code=code,cap=cap,seed=seed,**timing,search_evals=sum(int(a)+int(b) for a,b in parts),
                      explicit_evals=sum(map(int,refined))))
    write_csv(o/"cost.csv",c)
    (o/"binary.json").write_text(json.dumps(dict(sha256=hashlib.sha256(binary.read_bytes()).hexdigest()),indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("sugar",type=Path)
    p.add_argument("out",type=Path)
    p.add_argument("binary",type=Path)
    a=p.parse_args()
    prepare(a.sugar.resolve(),a.out.resolve())
    dock(a.sugar.resolve(),a.out.resolve(),a.binary.resolve())
