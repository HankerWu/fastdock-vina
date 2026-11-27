"""Diagnostic: swap rigid placement and internal geometry to test two inaccessible native basins."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile
import time

import numpy as np
from scipy.spatial.distance import cdist, pdist

STAGE = Path(__file__).resolve().parents[2]/"directional_polar_validation/first_stage"


def coords(text):
    return {int(l[6:11]): np.array([float(l[i:i+8]) for i in (30,38,46)])
            for l in text.splitlines() if l.startswith(("ATOM  ","HETATM"))}


def stamp(template, xyz):
    lines = []
    visited = set()
    for l in template.splitlines():
        if l.startswith(("ATOM  ","HETATM")):
            serial = int(l[6:11])
            l = l[:30] + "".join(f"{x:8.3f}" for x in xyz[serial]) + l[54:]
            visited.add(serial)
        lines.append(l)
    assert visited==set(xyz)
    return "\n".join(lines)+"\n"


def rmsd(x, ref):
    return np.sqrt(((x-ref)**2).sum(axis=-1).mean(axis=-1)).min()


def transform(x, target):
    a,b=x.mean(0),target.mean(0)
    u,_,v=np.linalg.svd((x-a).T@(target-b))
    u[:,-1] *= np.sign(np.linalg.det(u@v))
    rot=u@v
    assert abs(np.linalg.det(rot)-1.)<1e-12
    return rot,b-a@rot


def energy(log):
    return sum(float(re.search(s+r"\s*:\s*([-\d.]+)", log)[1])
               for s in ("Final Intermolecular Energy", "Final Total Internal Energy"))


def reindexings(mappings):
    # These cognate depositions have one ligand altloc; fixed observed positions
    # are permuted by chemical graph automorphisms from the validated preparation.
    perms=[]
    for mapping in mappings:
        distances=cdist(mapping,mappings[0])
        indices=distances.argmin(axis=1)
        assert len(set(indices))==len(indices) and distances.min(axis=1).max()<1e-10
        perms.append(indices)
    return perms


def write(path, rows):
    with path.open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=rows[0],lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)


def main(out, binary, probe):
    source=out/"source"
    source.mkdir(parents=True)
    with tarfile.open(STAGE/"raw.tar.xz") as f:
        f.extractall(source,filter="data")
    samples=out/"runs"
    samples.mkdir()
    results=[]
    cache={}
    for pdb in ("2YIY","4I1R"):
        site=source/pdb
        original=site/"ligand.pdbqt"
        template=original.read_text()
        receptor=site/"receptor.pdbqt"
        info=json.loads((site/"preparation.json").read_text())
        serial=info["serial"]
        ref=np.array(info["mapped"])
        perms=reindexings(ref)
        projected=coords((source/"crystal_diagnostic"/f"{pdb}-{pdb}-best_projected-input.pdbqt").read_text())
        px=np.array([projected[i] for i in serial])
        box=[]
        for axis,center in zip("xyz",info["center"]):
            box+=["--center_"+axis,str(center),"--size_"+axis,"48"]
        def evaluate(site_xyz):
            start=stamp(template,site_xyz)
            key=pdb,hashlib.sha256(start.encode()).hexdigest()
            if key not in cache:
                base=samples/(pdb+"-"+key[1][:16])
                input_path=base.with_suffix(".pdbqt")
                final=base.with_suffix(".out.pdbqt")
                input_path.write_text(start)
                cli=[str(binary),"--receptor",str(receptor),"--cpu","1","--seed","43103"]
                before=time.perf_counter()
                log=subprocess.check_output(cli+["--ligand",str(input_path),"--local_only","--out",str(final),*box],text=True)
                local_s=time.perf_counter()-before
                base.with_suffix(".log").write_text(log)
                checks=[]
                for path in (input_path,final):
                    log2=subprocess.check_output(cli+["--ligand",str(original),"--pose_coordinates",str(path),"--score_only"],text=True)
                    path.with_suffix(".score.log").write_text(log2)
                    checks.append(energy(log2))
                topo=subprocess.check_output([str(probe)],input=f"{original}\t{input_path}\n{original}\t{final}\n",text=True).splitlines()
                values=[line.split() for line in topo]
                assert all(x[:5]==["0"]*5 for x in values), (pdb,values)
                get_heavy=lambda content: np.array([coords(content)[j] for j in serial])
                cache[key]=dict(artifact=base.name, input_max_abs_offset_A=float(np.abs(get_heavy(start)-info["center"]).max()),input_rmsd=rmsd(get_heavy(start),ref),
                                final_rmsd=rmsd(get_heavy(final.read_text()),ref),input_energy=checks[0],final_energy=checks[1],
                                local_log_energy=energy(log),local_wall_s=local_s,max_bond_delta_A=max(float(x[-1]) for x in values),
                                max_pairdist_rounding_A=np.max(abs(pdist(get_heavy(start))-pdist(np.array([site_xyz[i] for i in serial])))))
            return cache[key]
        for path in sorted((source/"runs").glob(pdb+"_*_pro.pdbqt")):
            _,cap,seed,_=path.stem.split("_")
            blind=[coords(b) for b in path.read_text().split("ENDMDL") if b.strip()]
            rms=[rmsd(np.array([m[i] for i in serial]),ref) for m in blind]
            best=int(np.argmin(rms))
            for idx in dict.fromkeys([0,best]):
                pose=blind[idx]
                bx=np.array([pose[i] for i in serial])
                fits=[]
                for k,perm in enumerate(perms):
                    rot,translation=transform(bx,px[perm])
                    d=np.sqrt(((bx@rot+translation-px[perm])**2).sum(axis=1).mean())
                    fits.append((d,k,rot,translation))
                d,k,rot,translation=min(fits,key=lambda item:(item[0],item[1]))
                interventions={"blind":pose,"blind_shape_projected_frame":{i:v@rot+translation for i,v in pose.items()},
                               "projected_shape_blind_frame":{i:(v-translation)@rot.T for i,v in projected.items()},
                               "projected":projected}
                for arm,xyz in interventions.items():
                    rows=evaluate(xyz)
                    results.append(dict(pdb=pdb,cap=cap,seed=seed,rank=idx+1,
                                        selection="top_and_closest" if best==0 else ("top" if idx==0 else "closest"),
                                        symmetry=k,crossfit_internal_rmsd=d,arm=arm,**rows))
    write(out/"factorial.csv",results)
    meta=dict(vina_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),local_refinements=len(cache),rows=len(results),
              total_local_wall_s=sum(s["local_wall_s"] for s in cache.values()),
              max_rounding_pair_distance_change=max(s["max_pairdist_rounding_A"] for s in cache.values()),
              largest_local_vs_fixedtopo_output_energy_abs=max(abs(s["final_energy"]-s["local_log_energy"]) for s in cache.values()))
    (out/"audit.json").write_text(json.dumps(meta,indent=2)+"\n")
    print(meta)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ("work","vina","probe"):
        p.add_argument(flag,type=Path)
    a=p.parse_args()
    main(a.work.resolve(),a.vina.resolve(),a.probe.resolve())
