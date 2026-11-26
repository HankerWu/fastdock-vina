"""Check the observed alternative BRDT pocket copies on the same Kd pair."""
import argparse
import csv
import hashlib
import json
import math
import shutil
import statistics
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "geometry_metric"))
import panel
sys.path.insert(0, str(HERE.parents[1]))
from pose_run import align_with_numbering, load_meta, measure, write_csv

PDB = {"7MRH":"ZMJ","7MRG":"ZMM"}
SEEDS = (95701,95713,95717)


def prepare_b(inputs, work, code, ligand):
    original = inputs/code/(code+".pdb")
    directory = work/code
    directory.mkdir(parents=True, exist_ok=True)
    # The existing preparation chooses the first lexical instance. Interchange
    # only A and B labels, retaining every deposited atom for its exclusion tests.
    source = original.read_text().splitlines()
    swap = {"A":"B","B":"A"}
    changed = [l[:21]+swap[l[21]]+l[22:] if l.startswith(("ATOM  ","HETATM","TER   ")) and len(l)>21
               and l[21] in swap else l for l in source]
    (directory/(code+".pdb")).write_text("\n".join(changed)+"\n")
    shutil.copyfile(inputs/code/(ligand+".cif"), directory/(ligand+".cif"))
    report = {"pdb": code, "auth_chain":"B", "original_pdb_sha256":hashlib.sha256(original.read_bytes()).hexdigest()}
    try:
        meta = panel.prepare(work,code,ligand,"BRDT")
        if meta["chain"] != "A" or (directory/"ligand.pdbqt").read_bytes() != (inputs/code/"ligand.pdbqt").read_bytes():
            raise ValueError("Expected original chain B and the identical graph-prepared ligand")
        report.update(status="ready", nearest_omitted_nonwater_A=meta["nearest_omitted_nonwater_A"],
                      nearest_repaired_heavy_A=meta["closest_repaired_ligand_distance"])
    except (ValueError, KeyError, subprocess.CalledProcessError) as exc:
        report.update(status="excluded", reason=(str(exc) if not isinstance(exc,subprocess.CalledProcessError)
                     else "Meeko receptor template matching failed. See prep.log."))
        source_meta = directory/"source.json"
        if source_meta.is_file():
            report["nearest_omitted_nonwater_A"] = json.loads(source_meta.read_text())["nearest_omitted_nonwater_A"]
    return report


def ligand_only(folder, code, ligand):
    # A receptor-template failure does not invalidate its observed ligand or backbone alignment.
    d=folder/code
    lines=(d/(code+".pdb")).read_text().splitlines()
    refs=panel.base.dep_records(lines,ligand)
    meta=panel.base.ligand(d/(ligand+".cif"),refs,d/"ligand.pdbqt")
    meta.update(pdb=code,target="BRDT",chain="A",center=np.mean(list(refs[0].values()),axis=0).tolist())
    meta["ca"]={l[22:27].strip()+":"+l[17:20]:coordinates(l)
                for l in lines if l.startswith("ATOM") and l[21]=="A"
                and l[12:16].strip()=="CA" and l[16] in (" ","A")}
    (d/"query_only.json").write_text(json.dumps(meta,indent=2,default=lambda x:x.tolist())+"\n")
    return meta


def coordinates(line):
    return [float(line[i:i+8]) for i in (30,38,46)]


def lys_nz(lines, chain):
    hits = [coordinates(l) for l in lines if l.startswith("ATOM") and l[21] == chain and l[17:20]=="LYS"
            and int(l[22:26])==60 and l[12:16].strip()=="NZ" and l[16] in (" ","A")]
    if len(hits)!=1:
        raise ValueError("Ambiguous Lys60 NZ")
    return hits[0]


def contacts(raw,ranked,lys):
    by_rank = {}
    for i,model in enumerate(raw.split("ENDMDL")):
        oxygen = [coordinates(l) for l in model.splitlines() if l.startswith(("ATOM  ","HETATM"))
                  and int(l[6:11]) in (4,5) and l[77:].strip() == "OA"]
        if oxygen:
            assert len(oxygen)==2
            by_rank[i+1] = min(math.dist(lys,x) for x in oxygen)
    result = {}
    for label,r in [("first",ranked[0]),("nearest",min(ranked,key=lambda r:float(r["rmsd_A"])))]:
        rank=int(r["rank"])
        result.update({label+"_rank":rank,label+"_rmsd_A":float(r["rmsd_A"]),label+"_lys60_A":by_rank[rank]})
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("original",type=Path,help="Extracted water transfer inputs, dry outputs, modes and cells")
    p.add_argument("out",type=Path)
    p.add_argument("binary",type=Path)
    p.add_argument("--jobs",type=int,default=1)
    a=p.parse_args()
    inputs=a.original.resolve()/"inputs"
    work=a.out.resolve(); b=work/"B"
    manifest=[prepare_b(inputs,b,code,ligand) for code,ligand in PDB.items()]
    (work/"manifest_B.json").write_text(json.dumps(manifest,indent=2)+"\n")
    ready={c["pdb"]:load_meta(b,c["pdb"]) for c in manifest if c["status"]=="ready"}
    rows=[]; modes=[]
    for name,table in [("cells.csv",rows),("modes.csv",modes)]:
        with (a.original/name).open() as stream:
            table.extend(dict(r,chain="A") for r in csv.DictReader(stream) if r["arm"]=="dry")
    jobs=[]
    out=b/"poses";out.mkdir(exist_ok=True)
    queries=dict(ready)
    for code,ligand in PDB.items():
        if code not in queries:
            queries[code]=ligand_only(b,code,ligand)
        if (b/code/"ligand.pdbqt").read_bytes()!=(inputs/code/"ligand.pdbqt").read_bytes():
            raise ValueError("Prepared input differs between chains")
    for q,qm in queries.items():
        for r,rm in ready.items():
            rot,tr,rms,n,offset=np.eye(3),np.zeros(3),0.,0,0
            if q!=r:rot,tr,rms,n,offset=align_with_numbering(qm,rm)
            ref=np.asarray(qm["mapped"])@rot+tr
            for seed in SEEDS:
                jobs.append((b,out,a.binary.resolve(),qm,rm,ref,rms,n,offset,seed,"dry",120000))
    with ThreadPoolExecutor(a.jobs) as pool:
        for row,ms in pool.map(measure,jobs):
            rows.append(dict(row,chain="B")); modes.extend(dict(m,chain="B") for m in ms)
    write_csv(work/"cells.csv",rows)
    write_csv(work/"modes.csv",modes)
    contact_rows=[]
    for row in rows:
        if row["query"]!="7MRH":continue
        r,s,c=row["receiver"],str(row["seed"]),row["chain"]
        d=(a.original/"dry") if c=="A" else b
        raw=(d/"poses"/f"7MRH-{r}-{s}-dry.pdbqt").read_text()
        ranked=[m for m in modes if m["query"]=="7MRH" and m["receiver"]==r and str(m["seed"])==s and m["chain"]==c]
        lys=lys_nz((d/r/"receptor.pdbqt").read_text().splitlines(),"A")
        contact_rows.append(dict(receiver=r,chain=c,seed=s,**contacts(raw,ranked,lys)))
    write_csv(work/"contacts.csv",contact_rows)
    order=[]
    for c in ("A","B"):
        for r in PDB:
            selected=[x for x in rows if x["chain"]==c and x["receiver"]==r]
            if {x["query"] for x in selected}==set(PDB):
                med={q:statistics.median(float(x["score"]) for x in selected if x["query"]==q) for q in PDB}
                order.append(dict(receiver=r,chain=c,ZMJ_77nM=med["7MRH"],ZMM_875nM=med["7MRG"],
                                  correct=med["7MRH"]<med["7MRG"]))
    write_csv(work/"orders.csv",order)


if __name__=="__main__":
    main()
