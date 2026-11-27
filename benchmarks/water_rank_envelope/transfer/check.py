"""One registered AChE cross-dock and its frozen water-model certificates."""
import argparse
import hashlib
import json
import lzma
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import tarfile
import time

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
sys.path.insert(0, str(BENCH / "paired_crossdocking"))
import run as chemistry
sys.path.insert(0, str(HERE.parent))
from study import Site, availability, cross_sites, models, pdb_lines, penalty, xyz
from envelope import assess
from mixtures import columns, game, TOL


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tick():
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return np.array([time.perf_counter(), time.process_time() + children.ru_utime + children.ru_stime])


def ca(rows):
    return {line[22:27] + line[17:20]: xyz(line) for line in rows
            if line.startswith("ATOM") and line[21] == "A" and line[12:16] == " CA "}


def coords(lines, order):
    records = {int(line[6:11]): xyz(line) for line in lines if line.startswith(("ATOM", "HETATM"))}
    return np.asarray([records[serial] for serial in order])


def main(work, out, executable, arm):
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    clock = tick()
    sources = BENCH / "receptor_amide_states/transfer/transfer_evidence.tar.xz"
    with tarfile.open(sources) as source:
        for member in ("raw/4EY6.pdb", "4EY6/receptor.pdbqt", "4EY6/preparation.json"):
            (work/Path(member).name).write_bytes(source.extractfile("novel/" + member).read())
    (work/"4EY5.pdb").write_bytes(lzma.decompress((HERE/"4EY5.pdb.xz").read_bytes()))
    rec = json.loads((work/"preparation.json").read_text())
    assert rec["receptor_status"] == "ready" and rec["metals"] == []
    assert rec["closest_omitted_heavy_atom"][0] >= 8 and not rec["missing_atoms"]
    assert sha(work/"4EY6.pdb") == rec["pdb_sha256"]
    a, b = [pdb_lines(work/f"{p}.pdb") for p in ("4EY5", "4EY6")]
    crystal = chemistry.dep_records(a,"HUP")
    # This is the only input consumed by the docking binary.
    query = chemistry.ligand(HERE/"HUP.cif", [], work/"HUP.pdbqt")
    # Build the evaluation correspondence into a distinct file and compare bytes.
    prepared = chemistry.ligand(HERE/"HUP.cif", crystal, work/"mapped_reference_check.pdbqt")
    graph_only_digest = sha(work/"HUP.pdbqt")
    if sha(work/"mapped_reference_check.pdbqt") != graph_only_digest:
        raise ValueError("Supplying evaluation reference changed the input coordinates")
    assert query["mapped"] == [] and query["serial"] == prepared["serial"]
    qcenter = np.mean(list(crystal[0].values()),axis=0)
    rcenter = np.array(rec["center"])
    center_from_gnt = np.mean(list(chemistry.dep_records(b,"GNT")[0].values()),axis=0)
    assert np.allclose(rcenter, center_from_gnt, atol=1e-10, rtol=0)
    rot, trans, ca_rmsd, nca = chemistry.alignment(dict(ca=ca(a),center=qcenter),
                                                    dict(ca=ca(b),center=rcenter))
    reference = np.asarray(prepared["mapped"]) @ rot + trans
    molecule = Chem.MolFromSmiles(prepared["canonical_smiles"])
    assert reference.shape[1] == molecule.GetNumHeavyAtoms()
    costs = {"load_and_prepare_wall_cpu_s": (tick()-clock).tolist()}
    clock = tick()
    command = ["vina", "--receptor", "receptor.pdbqt", "--ligand", "HUP.pdbqt",
               "--cpu","1", "--seed","1103501","--exhaustiveness","3","--max_evals","24000",
               "--num_modes","9","--energy_range","100","--search_counters","--out","dry.pdbqt"]
    for axis, center in zip("xyz", rcenter, strict=True):
        command += ["--center_"+axis, str(center), "--size_"+axis, "24"]
    stored = HERE if arm == "pro24k" else HERE / arm
    if arm == "pro48k":
        command[command.index("--max_evals")+1] = "48000"
    if arm == "official24k":
        command.remove("--search_counters")
        if executable and sha(executable) != "f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644":
            raise ValueError("Incorrect official 1.2.7 executable digest")
    if executable:
        run = subprocess.run([str(executable.resolve()), *command[1:]], cwd=work,
                             capture_output=True, text=True)
        (out/"dry.log").write_text(run.stdout+run.stderr)
        run.check_returncode()
        shutil.copy(work/"dry.pdbqt",out/"dry.pdbqt")
        shutil.copy(work/"HUP.pdbqt",out/"HUP.pdbqt")
        costs["dock_wall_cpu_s"] = (tick()-clock).tolist()
    else:
        if out.resolve() == stored.resolve():
            raise ValueError("Use a fresh output directory for offline verification")
        shutil.copy(stored/"dry.pdbqt",work/"dry.pdbqt")
        assert sha(work/"HUP.pdbqt") == sha(HERE/"HUP.pdbqt")
    clock = tick()
    trees = list(models(work/"dry.pdbqt"))
    blocks = (work/"dry.pdbqt").read_text().strip().split("ENDMDL")[:-1]
    d = np.array([np.min(np.sqrt(np.square(coords(block.splitlines(),prepared["serial"])-reference)
                                .sum(axis=-1).mean(axis=1))) for block in blocks])
    scores = np.array([s for _,s,_ in trees])
    assert len(d)==len(trees) and len(trees)>0
    protein = cKDTree([xyz(l) for l in b if l.startswith("ATOM") and l[21]=="A"])
    records = cross_sites(b,protein,rcenter)
    site_rows = [dict(site=i+1, residue=code, x=c[0], y=c[1], z=c[2], B=B)
                 for i,(code,c,B) in enumerate(records)]
    pd.DataFrame(site_rows).to_csv(out/"sites.csv",index=False)
    factor_rows, rows, summary = [], [], []
    for order in (11,17):
        sites = [Site.prepare(p,B,protein,order=order) for _,p,B in records]
        avail = np.array([[availability(site,t) for site in sites] for _,_,t in trees])
        water = scores + penalty(avail).sum(axis=1)
        details, lower, upper, mask, winners = assess(scores,avail,d)
        records_for_order = []
        # Sites identically equal to one for every pose cannot affect the ranking.
        active = np.ptp(avail,axis=0)>0
        logs = None
        for i, (rival, margin) in enumerate(details):
            status = "single_rival" if rival is not None else "open"
            if status=="open" and lower < upper-1e-10:
                if logs is None:
                    if active.sum() > 20:
                        raise ValueError("Exact joint LP exceeds registered enumeration capacity")
                    logs=columns(scores,avail[:,active])
                best, worst = game(logs,i,+1), game(logs,i,-1)
                if best["upper"] < -TOL:
                    status = "mixed_rival"
                elif worst["lower"] > TOL:
                    status = "realizable"
            records_for_order.append(dict(order=order, mode=i+1, dry=scores[i], hydrated=water[i],
                                          rmsd=d[i], status=status, rival=None if rival is None else rival+1,
                                          margin=None if not np.isfinite(margin) else margin))
            for j,value in enumerate(avail[i],1):
                factor_rows.append(dict(order=order,mode=i+1,site=j,available=value))
        for row in records_for_order:
            if row["status"]=="realizable":
                upper=min(upper,row["rmsd"])
        lower=min(row["rmsd"] for row in records_for_order if row["status"] not in ("single_rival","mixed_rival"))
        rows.extend(records_for_order)
        dry_rank,wet_rank=np.argsort(scores,kind="stable"),np.argsort(water,kind="stable")
        io=int(np.argmin(d))
        summary.append(dict(order=order, sites=len(sites), active_sites=int(active.sum()), modes=len(d),
                            dry_mode=int(dry_rank[0])+1, water_mode=int(wet_rank[0])+1, oracle_mode=io+1,
                            dry=d[dry_rank[0]], water=d[wet_rank[0]], oracle=d[io],
                            dry_best5=min(d[dry_rank[:5]]), water_best5=min(d[wet_rank[:5]]),
                            dry_regret=d[dry_rank[0]]-d[io],water_regret=d[wet_rank[0]]-d[io],
                            lower=lower,upper=upper, excluded=sum(r["status"] in ("single_rival","mixed_rival")
                                                               for r in records_for_order),
                            witness_initial_sites=[j+1 for j,yes in enumerate(mask) if yes],
                            witness_initial_modes=[int(i)+1 for i in winners]))
    for name, data in (("modes",rows),("factors",factor_rows),("summary",summary)):
        pd.DataFrame(data).to_csv(out/(name+".csv"),index=False)
    costs["site_evaluation_wall_cpu_s"]=(tick()-clock).tolist()
    info=dict(arm=arm,command=command, binary_sha256=sha(executable) if executable else None,
              source_archive_sha256=sha(sources), source_urls=["https://files.rcsb.org/download/4EY5.pdb",
                "https://files.rcsb.org/ligands/download/HUP.cif","https://www.rcsb.org/structure/4EY6"],
              hashes={p.name:sha(p) for p in (work/"4EY5.pdb",HERE/"HUP.cif",work/"4EY6.pdb",
                work/"receptor.pdbqt",work/"preparation.json",work/"HUP.pdbqt",work/"dry.pdbqt")},
              graph_only_pdbqt_sha256=graph_only_digest, matched_reference_pdbqt_sha256=sha(work/"mapped_reference_check.pdbqt"),
              box_center_from_receiver_gnt=center_from_gnt.tolist(), smiles=prepared["canonical_smiles"],formal_charge=Chem.GetFormalCharge(molecule),
              atoms=molecule.GetNumHeavyAtoms(), mappings=len(reference), receptor_alignment_A=ca_rmsd,
              nca=nca, local_metal=rec["metals"], nearest_omitted_A=rec["closest_omitted_heavy_atom"][0],
              costs=costs)
    if not executable:
        archive = json.loads((stored/"manifest.json").read_text())
        if info["hashes"] != archive["hashes"] or info["source_archive_sha256"] != archive["source_archive_sha256"]:
            raise ValueError("Archived input or pose bytes changed")
    (out/"manifest.json").write_text(json.dumps(info,indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("out",type=Path)
    p.add_argument("--vina",type=Path)
    p.add_argument("--arm",choices=("pro24k","pro48k","official24k"),default="pro24k")
    a=p.parse_args()
    main(a.work,a.out,a.vina,a.arm)
