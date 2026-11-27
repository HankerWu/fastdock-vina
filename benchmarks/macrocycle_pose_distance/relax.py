"""Relax common-topology original and rigidly placed poses with the native local optimizer."""
import argparse
import csv
import hashlib
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "macrocycle_proposal" / "external"))
from external_score import poses, reference, write


def bare(source, dest):
    dest.write_text("\n".join(l for l in source.read_text().splitlines()
                              if not l.startswith(("MODEL ", "ENDMDL")))+"\n")


def frozen_score(engine, folder, box, pose, log):
    cmd = [engine, "--score_only", "--pose_coordinates", pose, "--ligand", folder/"ligand.pdbqt",
           "--receptor", folder/"receptor.pdbqt", "--config", box, "--cpu", "1", "--seed", "458921",
           "--unbound_energy", "0"]
    output = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
    log.write_text(output)
    terms = re.findall(r":\s*([-+.\deE]+)\s+\(kcal/mol\)", output)
    if len(terms) != 10:
        raise ValueError(f"Incomplete score for {pose}")
    return {"energy":float(terms[0]),"inter":float(terms[1]),"intra":float(terms[4])}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs","results","placement","engine","probe","out"):
        p.add_argument(name,type=Path)
    p.add_argument("--observer",type=Path)
    a=p.parse_args()
    a.out.mkdir(exist_ok=True,parents=True)
    by_structure,path_rows={},[]
    keys=sorted(path.name.removesuffix("_native.pdbqt") for path in a.placement.glob("*_native.pdbqt"))
    for key in keys:
        pdb=key[:4]
        folder=a.inputs/pdb
        if pdb not in by_structure:
            by_structure[pdb] = reference(a.inputs,pdb)
        serial, refs, order, sdf, ends, _, tree, bond = by_structure[pdb]
        box=a.results/(pdb+".box")
        target_sig=subprocess.check_output([a.probe.resolve(),folder/"ligand.pdbqt"])
        for arm in ("native","rigid"):
            stem=key+"_"+arm
            start=a.out/(stem+".start.pdbqt")
            bare(a.placement/(stem+".pdbqt"),start)
            start_sig=subprocess.check_output([a.probe.resolve(),start])
            if target_sig!=start_sig:
                raise ValueError("Changed native topology: " + stem)
            dest=a.out/(stem+".local.pdbqt")
            args=["--local_only","--ligand",start,"--receptor",folder/"receptor.pdbqt","--config",box,
                  "--cpu","1","--seed","458921","--verbosity","2","--out",dest]
            output=subprocess.check_output([a.engine.resolve(),*args],text=True,stderr=subprocess.STDOUT)
            (a.out/(stem+".local.log")).write_text(output)
            steps=int(re.search(r"Number of local optimization steps: (\d+)",output)[1])
            calls=""
            if a.observer:
                observed=a.out/(stem+".observed.pdbqt")
                trace=subprocess.check_output([a.observer.resolve(),*args[:-1],observed],text=True,stderr=subprocess.STDOUT)
                assert observed.read_bytes()==dest.read_bytes()
                observed.unlink()
                (a.out/(stem+".trace")).write_text(trace)
                calls=int(re.search(r"NATIVE_LOCAL_EVALUATIONS (\d+)",trace)[1])
            for stage,pose in (("before",start),("after",dest)):
                scores=frozen_score(a.engine.resolve(),folder,box,pose,a.out/(stem+"."+stage+".score"))
                xyz=next(poses(pose))
                chem=np.array([xyz[i][1] for i in serial])
                coords=np.array([xyz[i][1] for i in order])
                rms=np.sqrt(((chem-refs)**2).sum(axis=-1).mean(axis=-1)).min()
                rms2=np.sqrt(((coords-sdf)**2).sum(axis=-1).mean(axis=-1)).min()
                if abs(rms-rms2)>0.0015:
                    raise ValueError(f"Reference RMSD disagreement for {stem}/{stage}")
                nearest=tree.query(chem)[0]
                gap=np.linalg.norm(xyz[ends[0]][1]-xyz[ends[1]][1])
                path_rows.append(dict(key=key,arm=arm,stage=stage,steps_limit=steps,
                    actual_gradient_evals=calls if stage=="after" else "",rmsd_A=rms,sdf_rmsd_A=rms2,
                    clash_atoms=sum(nearest<2),nearest_A=nearest.min(),cg_A=gap,native_cg_A=bond,
                    centroid_error_A=np.linalg.norm(chem.mean(axis=0)-refs[0].mean(axis=0)),**scores,
                    topology_sha256=hashlib.sha256(start_sig).hexdigest()))
    write(a.out/"relaxation.csv",path_rows)
