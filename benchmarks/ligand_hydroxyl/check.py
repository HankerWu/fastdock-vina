"""Rotate only observed ligand polar H; independently evaluate native and torsion energies."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation

from model import (BENCH, Receptor, atoms, hydrogenated_positions, ring, score, topology,
                   torsional_prior, unit)
from evaluate import archive_cases, blocks, write  # archived common input iterator


def rotate(block, sites, degrees):
    xyz = atoms(block)
    changes = {}
    for site in sites.values():
        oxygen, carbon, hydrogen = (xyz[i].xyz for i in site[:3])
        axis = unit(oxygen-carbon)
        h = oxygen + Rotation.from_rotvec(np.deg2rad(degrees)*axis).apply(hydrogen-oxygen)
        changes[site.hydrogen] = "".join(f"{v:8.3f}" for v in h)
    rows = []
    for line in block.splitlines():
        if line.startswith(("ATOM  ","HETATM")) and int(line[6:11]) in changes:
            line = line[:30]+changes[int(line[6:11])]+line[54:]
        rows.append(line)
    return "\n".join(rows)+"\n"


def forcefield_scan(mol,serial,xyz,site):
    all_xyz = hydrogenated_positions(mol,serial,xyz)
    full = Chem.AddHs(mol)
    conf = Chem.Conformer(full.GetNumAtoms())
    for i,v in enumerate(all_xyz):
        conf.SetAtomPosition(i,tuple(v))
    full.AddConformer(conf)
    prop = AllChem.MMFFGetMoleculeProperties(full,mmffVariant="MMFF94s")
    for term in ("Bond","Angle","StretchBend","Oop","VdW","Ele"):
        getattr(prop, "SetMMFF"+term+"Term")(False)
    field = AllChem.MMFFGetMoleculeForceField(full,prop)
    oi = next(i for i,s in serial.items() if s==site.oxygen)
    h = next(a.GetIdx() for a in full.GetAtomWithIdx(oi).GetNeighbors() if a.GetAtomicNum()==1)
    nodes,around,axis = ring(site,xyz,72)
    expected = torsional_prior(site,around,axis,xyz,serial,all_xyz)
    observed = []
    for v in nodes:
        all_xyz[h] = v
        observed.append(field.CalcEnergy(all_xyz.ravel().tolist()))
    return float(np.max(np.abs((np.array(observed)-observed[0])-(expected-expected[0]))))


def verify(binary, out, revision):
    out.mkdir(parents=True,exist_ok=True)
    rows, scans = [], []
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        for base,sites,pool,_ in archive_cases(root):
            key = (base["panel"],base["design"],base["query"],base["seed"])
            controls = {("esr2","12000","5TOA","52711"),("esr2","60000","1QKM","52711"),
                        ("cognate","matched","7QK0","85237"),("reciprocal","fixed","1RBP","85237")}
            if key not in controls:
                continue
            receptor = sites/base["receiver"]/"receptor.pdbqt"
            ligand = sites/base["query"]/"ligand.pdbqt"
            rec = Receptor(receptor.read_text())
            query = ligand.read_text()
            mol,serial,oh = topology(query)
            for rank,block in enumerate(blocks(pool)[:2],1):
                xyz = atoms(block)
                for site in oh.values():
                    err = forcefield_scan(mol,serial,xyz,site)
                    if err>1.e-8:
                        raise AssertionError((key,rank,site.oxygen,err))
                    scans.append(dict(**base,rank=rank,oxygen=site.oxygen,error_kcal_mol=err))
                native_ref = None
                for shift in (0,180,37):
                    rotated = rotate(block,oh,shift)
                    coords = atoms(rotated)
                    if any(np.any(a.xyz!=coords[i].xyz) for i,a in xyz.items() if a.kind not in ("HD","H")):
                        raise ValueError("Heavy atom moved")
                    pose = root/"counterfactual.pdbqt"
                    pose.write_text(rotated+"ENDMDL\n")
                    cli = [str(binary),"--cpu","1","--seed","43103","--score_only","--receptor",str(receptor),
                           "--ligand",str(ligand),"--pose_coordinates",str(pose),"--net_polar","native","--verbosity","0"]
                    log = subprocess.check_output(cli,text=True)
                    native = next(l for l in log.splitlines() if l.startswith("Net polar native:")).split()[3:]
                    if native_ref is None:
                        native_ref = native
                    if native != native_ref:
                        raise AssertionError((key,native_ref,native))
                    terms = score(rec,query,rotated)
                    rows.append(dict(**base,rank=rank,shift=shift,native_inter=native[0],native_intra=native[1],
                                     normalized_zero_reference=native[2],raw_fixed=terms["fixed"],
                                     raw_bonded=terms["bonded"],raw_uniform=terms["uniform"]))
    if len(rows)!=24:
        raise AssertionError(len(rows))
    write(out/"hydrogen_replay.csv",rows)
    write(out/"mmff_scans.csv",scans)
    (out/"binary.json").write_text(json.dumps(dict(reported_source_revision=revision,
                                      sha256=hashlib.sha256(binary.read_bytes()).hexdigest()),indent=2)+"\n")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--revision", help="Source revision used to build this executable")
    args=parser.parse_args()
    verify(args.binary.resolve(),args.output,args.revision)
