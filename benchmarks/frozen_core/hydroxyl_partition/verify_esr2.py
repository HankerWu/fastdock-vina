"""Redock the same two preparations with upstream and Pro, and independently measure every pose."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import resource
import subprocess
import tarfile
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def xyz(line):
    return [float(line[i:i+8]) for i in (30, 38, 46)]


def audit_raw(work, native, pro):
    here = Path(__file__).resolve().parent
    with tarfile.open(here/"esr2.tar.xz") as archive:
        archive.extractall(work,filter="data")
    out = work/"control"
    out.mkdir()
    rows, poses, environment = [], [], []
    for code, ligand in (("1QKM","GEN"), ("5TOA","EST")):
        site=work/code
        meta=json.loads((site/"preparation.json").read_text())
        deposition=(work/"raw"/(code+".pdb")).read_text().splitlines()
        waters=[xyz(l) for l in deposition if l.startswith("HETATM") and l[17:20]=="HOH" and l[76:78].strip()=="O"]
        ligand_xyz=[xyz(l) for l in deposition if l.startswith("HETATM") and l[17:20]==ligand and l[21]=="A" and l[76:78].strip()!="H"]
        environment.append(dict(pdb=code,resolution_A=float(next(l for l in deposition if l.startswith("REMARK   2 RESOLUTION.")).split()[3]),
                                min_water_distance_A=min((math.dist(w,l) for w in waters for l in ligand_xyz),default=None),
                                min_other_omitted_heavy=meta["closest_omitted_heavy_atom"], min_repaired_A=meta["closest_repaired_ligand_distance"],
                                symmetry_mappings=len(meta["mapped"]), heavy_atoms=len(meta["serial"])))
        for record in sorted(site.glob("*-*-*.json")):
            args=json.loads(record.read_text())["cmd"]
            for method,binary in (("upstream",native),("pro",pro)):
                v=args.copy()
                v[0]=str(binary)
                output=out/(method+"-"+record.stem+".pdbqt")
                for flag in ("--receptor","--ligand"):
                    i=v.index(flag)+1
                    v[i]=str(site/v[i])
                v[v.index("--out")+1]=str(output)
                before=resource.getrusage(resource.RUSAGE_CHILDREN)
                start=time.perf_counter()
                with output.with_suffix(".log").open("w") as f:
                    subprocess.run(v,stdout=f,stderr=subprocess.STDOUT,check=True)
                duration=time.perf_counter()-start
                after=resource.getrusage(resource.RUSAGE_CHILDREN)
                byte_equal=output.read_bytes()==record.with_suffix(".pdbqt").read_bytes()
                if method=="pro" and not byte_equal:
                    raise ValueError("Rerun disagrees with archived Pro pool")
                cell=dict(pdb=code,seed=v[v.index("--seed")+1],cap=v[v.index("--max_evals")+1],method=method)
                rms=[]
                for model in (b for b in output.read_text().split("ENDMDL") if b.strip()):
                    atoms={int(l[6:11]):xyz(l) for l in model.splitlines() if l.startswith(("ATOM  ","HETATM"))}
                    pos=[atoms[i] for i in meta["serial"]]
                    rmsd=min(math.sqrt(math.fsum(math.dist(a,b)**2 for a,b in zip(pos,m,strict=True))/len(pos)) for m in meta["mapped"])
                    rms.append(rmsd)
                    poses.append(dict(**cell,rank=len(rms),rmsd_A=rmsd))
                rows.append(dict(**cell,poses=len(rms),first=rms[0],best5=min(rms[:5]),best=min(rms),equal_archived_pro=byte_equal,
                                 wall_s=duration,child_cpu_s=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                                 binary_sha256=digest(binary),pool_sha256=digest(output)))
    for path,data in (("esr2_controls.csv",rows),("esr2_control_modes.csv",poses),("esr2_environment.csv",environment)):
        with (out/path).open("w",newline="") as f:
            wr=csv.DictWriter(f,fieldnames=data[0])
            wr.writeheader()
            wr.writerows(data)
    with (here/"results/modes.csv").open() as stream:
        observed={(r["query"],r["seed"],r["cap"],int(r["rank"])):float(r["rmsd_A"])
                  for r in csv.DictReader(stream) if r["panel"]=="esr2"}
    checks=[abs(observed[r["pdb"],r["seed"],r["cap"],r["rank"]]-r["rmsd_A"]) for r in poses if r["method"]=="pro"]
    assert len(checks)==len(observed) and max(checks)<1e-12
    print("Verified Pro modes:", len(checks), "worst RMSD difference:", max(checks))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("upstream",type=Path)
    p.add_argument("pro",type=Path)
    a=p.parse_args()
    audit_raw(a.work.resolve(),a.upstream.resolve(),a.pro.resolve())
