"""Paired carboxylate state and upstream force-field experiment in RBP4 and FABP4."""
import argparse
import itertools
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"confirmation"))
from run import measure, reference, write_csv

SERIES = {"2NNQ": ("2NNQ", "3FR2"), "3FR2": ("2NNQ", "3FR2"), "3FMZ": ("1RBP", "3FMZ")}
SEEDS = (93133, 93179, 93229)


def atoms(text):
    return [line for line in text.splitlines() if line.startswith(("ATOM  ","HETATM"))]


def encode(mol):
    (setup,) = MoleculePreparation().prepare(mol)
    text, ok, error = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(error)
    return text


def indices(text, n):
    result = {}
    for line in text.splitlines():
        if line.startswith("REMARK INDEX MAP"):
            values = list(map(int, line.split()[3:]))
            result.update((i-1,j) for i,j in zip(values[::2],values[1::2]) if i<=n)
    return [result[i] for i in range(n)]


def ligand(entry, fixtures, out):
    parent = Chem.MolFromSmiles(entry["canonical_smiles"])
    (_, oxygen1, oxygen2), = parent.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]"))
    mol=Chem.AddHs(parent)
    params=AllChem.ETKDGv3()
    params.randomSeed=43103
    if AllChem.EmbedMolecule(mol,params) or AllChem.MMFFOptimizeMolecule(mol,maxIters=1000):
        raise ValueError("Failed to reproduce neutral geometry")
    neutral=encode(mol)
    if neutral != (fixtures/entry["pdb"]/"ligand.pdbqt").read_text():
        raise ValueError("Neutral fixture differs")
    anion=Chem.RWMol(mol)
    h, = [a.GetIdx() for a in anion.GetAtomWithIdx(oxygen2).GetNeighbors() if a.GetAtomicNum()==1]
    anion.GetAtomWithIdx(oxygen2).SetFormalCharge(-1)
    anion.RemoveAtom(h)
    Chem.SanitizeMol(anion)
    meta={}
    for state, text in (("HA",neutral),("A",encode(Chem.Mol(anion)))):
        (out/f"{entry['pdb']}-{state}.pdbqt").write_text(text)
        serial=indices(text,parent.GetNumAtoms())
        points={int(r[6:11]):r[30:54] for r in atoms(text)}
        xyz=[points[k] for k in serial]
        if state=="HA":
            baseline=xyz
        elif xyz != baseline:
            raise ValueError("Heavy atom coordinates changed")
        charge=sum(float(r[70:76]) for r in atoms(text))
        if abs(charge-(-1 if state=="A" else 0)) > 0.015:
            raise ValueError("Unexpected PDBQT net charge")
        meta[state]=dict(serial=serial,oxygen=[oxygen1,oxygen2],charge=charge,
                         torsdof=int(text.split("TORSDOF ")[1].split()[0]))
    return meta


def make_maps(program,fixtures,out,rec,atom_types):
    folder=out/rec["pdb"]
    folder.mkdir(exist_ok=True)
    text=(fixtures/rec["pdb"]/"receptor.pdbqt").read_text()
    (folder/"r.pdbqt").write_text(text)
    receptor_types=sorted({r.split()[-1] for r in atoms(text)})
    lines=["npts 64 64 64","gridfld r.maps.fld","spacing 0.375",
           "receptor_types "+" ".join(receptor_types),"ligand_types "+" ".join(atom_types),
           "receptor r.pdbqt","gridcenter "+" ".join(map(str,rec["center"])),"smooth 0.5"]
    lines += ["map r."+t+".map" for t in atom_types]+["elecmap r.e.map","dsolvmap r.d.map","dielectric -42"]
    (folder/"r.gpf").write_text("\n".join(lines)+"\n")
    subprocess.run([program,"-p","r.gpf","-l","r.glg"],cwd=folder,check=True,stdout=subprocess.DEVNULL)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ("fixtures","out","vina","autogrid"):
        ap.add_argument(name,type=Path)
    ap.add_argument("--jobs",type=int,default=2)
    a=ap.parse_args()
    fixtures,out,binary,autogrid=[p.resolve() for p in (a.fixtures,a.out,a.vina,a.autogrid)]
    out.mkdir(parents=True,exist_ok=True)
    entries={p.parent.name:json.loads(p.read_text()) for p in fixtures.glob("*/preparation.json")}
    states={q:ligand(entries[q],fixtures,out) for q in SERIES}
    (out/"states.json").write_text(json.dumps(states,indent=2)+"\n")
    types=sorted({r.split()[-1] for q,s in itertools.product(SERIES,("HA","A"))
                  for r in atoms((out/f"{q}-{s}.pdbqt").read_text())})
    for r in sorted(set(itertools.chain.from_iterable(SERIES.values()))):
        make_maps(autogrid,fixtures,out,entries[r],types)

    def dock(job):
        q,r,state,force,seed=job
        dest=out/f"{q}-{r}-{state}-{force}-{seed}.pdbqt"
        cmd=[str(binary),"--ligand",str(out/f"{q}-{state}.pdbqt"),"--cpu","1","--seed",str(seed),
             "--exhaustiveness","3","--max_evals","24000","--num_modes","9","--energy_range","100","--out",str(dest)]
        if force=="ad4":
            cmd.extend(("--scoring","ad4","--maps",str(out/r/"r")))
        else:
            cmd.extend(("--receptor",str(fixtures/r/"receptor.pdbqt"),"--spacing","0.375"))
            for axis,val in zip("xyz",entries[r]["center"]):
                cmd.extend(("--center_"+axis,str(val),"--size_"+axis,"24"))
        with dest.with_suffix(".log").open("w") as log:
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        info=states[q][state]
        expected=reference(entries[q],entries[r])
        oxy=info["oxygen"]
        permuted=expected.copy()
        permuted[:,oxy,:]=expected[:,oxy[::-1],:]
        plain=measure(dest,info["serial"],expected)
        equiv=measure(dest,info["serial"],np.concatenate([expected,permuted]))
        return [dict(target=entries[q]["target"],query=q,receiver=r,seed=seed,state=state,
                     force=force,rank=i,rmsd_fixed=fixed,rmsd_resonance=rms)
                for i,(fixed,rms) in enumerate(zip(plain,equiv),1)]
    jobs=[(q,r,s,f,k) for q,receivers in SERIES.items()
          for r,s,f,k in itertools.product(receivers,("HA","A"),("vina","ad4"),SEEDS)]
    with ThreadPoolExecutor(a.jobs) as runner:
        rows=[v for group in runner.map(dock,jobs) for v in group]
    write_csv(out/"modes.csv",rows)
    for q,s,f in itertools.product(SERIES,("HA","A"),("vina","ad4")):
        groups=[[v for v in rows if (v["query"],v["receiver"],v["seed"],v["state"],v["force"])==(q,r,k,s,f)]
                for r,k in itertools.product(SERIES[q],SEEDS)]
        print(q,s,f,*[sum(min(v["rmsd_resonance"] for v in g[:n])<=2 for g in groups) for n in (1,5,9)])


if __name__=="__main__":
    main()
