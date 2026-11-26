"""Frozen paired pilot on public crystal cases and source-stratified Ki series.

No query coordinates enter ligand preparation. Crystal-defined pocket boxes are
a known-site redock control, not a blind or cross-docking benchmark.
"""
import argparse
import csv
import itertools
import json
import math
import re
import subprocess
import time
from pathlib import Path

import numpy as np
import yaml
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from scipy.stats import spearmanr

POSE_IDS = """7WPW_F15 7NF3_4LU 8EX2_Q2Q 7N7H_CTP 7WJB_BGC 7UMW_NAD
7P1F_KFN 8AP0_PRP 5SB2_1K2 7B2C_TP7 7ED2_A3P 7VWF_K55 7NR8_UOE
6XBO_5MC 6TW5_9M2""".split()
AFFINITY = {"mcl1": "lig_43", "ptp1b": "lig_20670_2qbs", "tyk2": "lig_ejm_54"}
SEEDS = (20261126, 20261127, 20261128)
ARMS = {"native": (), "mixed": ("--mixed_vinardo_search",)}

def command(args, log, timeout=600):
    start = time.monotonic()
    p = subprocess.run([str(x) for x in args], capture_output=True, text=True, timeout=timeout)
    log.write_text(p.stdout + "\n" + p.stderr)
    if p.returncode:
        raise RuntimeError(f"exit {p.returncode}: {log}")
    return time.monotonic() - start

def read_sdf(path):
    mol = Chem.SDMolSupplier(str(path), removeHs=False)[0]
    if mol is None or mol.GetNumConformers() != 1:
        raise ValueError(f"invalid reference {path}")
    return mol

def prepare_ligand(smiles, path):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("invalid/multicomponent ligand")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 20261126
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError("ETKDG embedding failure")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    setups = MoleculePreparation().prepare(mol)
    if len(setups) != 1:
        raise ValueError("expected one ligand setup")
    text, ok, reason = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(reason)
    path.write_text(text)

def prepare_receptor(pdb, folder):
    prefix = folder / "receptor"
    out = prefix.with_suffix(".pdbqt")
    if not out.exists():
        # Fail closed: no deletion of unmatched residues/cofactors to rescue a case.
        command(["mk_prepare_receptor.py", "--read_pdb", pdb, "-o", prefix, "-p",
                 "--default_altloc", "A"], folder / "receptor.log")
    if not out.is_file() or not out.stat().st_size:
        raise ValueError("no rigid receptor output")
    return out

def coords(mol):
    return Chem.RemoveHs(mol).GetConformer().GetPositions()

def box_args(center):
    args = []
    for axis, value in zip("xyz", center):
        args += ["--center_" + axis, str(float(value)), "--size_" + axis, "25"]
    return args

def docking(binary, receptor, ligand, center, folder, seed, arm, exhaustiveness, cap):
    tag = f"{seed}_{arm}"
    out, log = folder / f"{tag}.pdbqt", folder / f"{tag}.log"
    flags = list(ARMS[arm])
    seconds = command([binary, "--receptor", receptor, "--ligand", ligand,
        *box_args(center), "--pose_selection", "late", "--exhaustiveness", exhaustiveness,
        "--max_evals", cap, "--refine_candidates", 40, "--num_modes", 10,
        "--energy_range", 100, "--cpu", 1, "--seed", seed, "--search_counters",
        *flags, "--out", out], log)
    text = log.read_text()
    scores = [float(x) for x in re.findall(r"REMARK VINA RESULT:\s+([-+0-9.eE]+)", out.read_text())]
    if not scores or not all(map(math.isfinite, scores)):
        raise ValueError("missing native score")
    if scores != sorted(scores):
        raise ValueError("native scores not ordered")
    if arm == "mixed" and "Mixed Vinardo search:" not in text:
        raise ValueError("mixed mode not reported")
    traj = re.findall(r"hunt_evals=(\d+) second_calls=\d+ second_evals=(\d+)", text)
    explicit = re.findall(r"explicit_evals=(\d+)", text)
    if len(traj) != exhaustiveness or not explicit:
        raise ValueError("missing evaluation counters")
    row = dict(score=scores[0], seconds=seconds, hunt_evals=sum(int(a) for a,b in traj),
               second_evals=sum(int(b) for a,b in traj), explicit_evals=int(explicit[-1]))
    return out, scores, row

def pose_metrics(out, scores, ref, pdb, folder):
    from posebusters import PoseBusters
    mols = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule.from_file(str(out), skip_typing=True))
    if len(mols) != 1 or mols[0] is None:
        raise ValueError("cannot restore ligand graph")
    poses = mols[0]
    if poses.GetNumConformers() != len(scores):
        raise ValueError("pose/score count mismatch")
    rmsds, valid, detail = [], [], []
    true = Chem.RemoveHs(ref)
    for rank, conf in enumerate(poses.GetConformers(), 1):
        # Explicitly preserve the receptor coordinate frame: CalcRMS does not align.
        p = Chem.Mol(poses)
        keep = Chem.Conformer(conf)
        p.RemoveAllConformers()
        p.AddConformer(keep, assignId=True)
        heavy = Chem.RemoveHs(p)
        if Chem.MolToSmiles(heavy) != Chem.MolToSmiles(true):
            raise ValueError("ligand graph/stereo changed")
        rms = float(rdMolAlign.CalcRMS(heavy, true, maxMatches=100000))
        # Independent graph-automorphism calculation (exact bond graph). The
        # relaxed RDKit terminal-group symmetry can only lower this RMSD.
        matches = true.GetSubstructMatches(heavy, uniquify=False, useChirality=True, maxMatches=100000)
        if not matches or len(matches) == 100000:
            raise ValueError("atom mapping absent or truncated")
        a, b = coords(heavy), coords(true)
        exact = min(float(np.sqrt(np.mean(np.sum((a-b[list(m)])**2, axis=1)))) for m in matches)
        check = float(rdMolAlign.CalcRMS(heavy, true, map=[list(enumerate(m)) for m in matches],
                                        symmetrizeConjugatedTerminalGroups=False))
        if abs(check-exact) > 1e-6 or rms > exact + 1e-6:
            raise ValueError("independent RMSD disagreement")
        sdf = folder / f"{out.stem}_pose{rank}.sdf"
        with Chem.SDWriter(str(sdf)) as writer:
            writer.write(p)
        df = PoseBusters(config="redock", max_workers=0).bust(str(sdf), mol_true=ref, mol_cond=str(pdb))
        if len(df) != 1:
            raise ValueError("unexpected PoseBusters row count")
        # Keep raw test verdicts. Separate RMSD from physical validity.
        tests = {str(k): (str(v).lower() == "true") for k,v in df.iloc[0].items()
                 if "rmsd" not in str(k).lower()}
        if len(tests) < 15:
            raise ValueError("unexpected PoseBusters test schema")
        ok = all(tests.values())
        rmsds.append(rms)
        valid.append(ok)
        detail.append(dict(rank=rank, rmsd=rms, exact_graph_rmsd=exact, pb_valid=ok, checks=tests))
        if rank == 5:
            break
    (folder / f"{out.stem}_metrics.json").write_text(json.dumps(detail, indent=2))
    return dict(top1_rmsd=rmsds[0], top5_rmsd=min(rmsds), top1_valid=valid[0],
                top1_success=rmsds[0] <= 2 and valid[0],
                top5_success=any(r <= 2 and v for r,v in zip(rmsds,valid)))

def pose_tasks(root, case_ids):
    paths = {p.stem.lower().removesuffix("_ligand"): p for p in
             root.rglob("*_ligand.sdf") if "posebusters_benchmark_set" in p.parts}
    for case in case_ids:
        ref = paths.get(case.lower())
        if ref is None:
            yield case, None, None
        else:
            pdb = ref.with_name(ref.name.replace("_ligand.sdf", "_protein.pdb"))
            yield case, ref, pdb

def run_pose(a):
    rows = []
    case_ids = POSE_IDS[int(a.subset)::3] if a.subset != 'all' else POSE_IDS
    for case, refpath, pdb in pose_tasks(a.data, case_ids):
        folder = a.output / case
        folder.mkdir(exist_ok=True, parents=True)
        common = dict(case=case)
        try:
            if refpath is None or not pdb.is_file():
                raise ValueError("frozen case missing from official archive")
            ref = read_sdf(refpath)
            receptor = prepare_receptor(pdb, folder)
            ligand = folder / "ligand.pdbqt"
            prepare_ligand(Chem.MolToSmiles(Chem.RemoveHs(ref)), ligand)
            center = coords(ref).mean(axis=0)
            for seed, arm in itertools.product(SEEDS, ARMS):
                row = dict(common, seed=seed, arm=arm, status="ok")
                try:
                    out, scores, values = docking(a.binary,receptor,ligand,center,folder,seed,arm,a.exhaustiveness,a.cap)
                    row.update(values)
                    row.update(pose_metrics(out,scores,ref,pdb,folder))
                except Exception as e:
                    row.update(status="failed", error=repr(e))
                rows.append(row)
        except Exception as e:
            rows.extend(dict(common,seed=s,arm=arm,status="prep_failed",error=repr(e))
                        for s,arm in itertools.product(SEEDS,ARMS))
        write_rows(a.output/"pose_rows.csv",rows)
    paired = []
    for case in case_ids:
        for seed in SEEDS:
            p = [r for r in rows if r["case"]==case and r["seed"]==seed]
            if len(p)==2 and all(r["status"]=="ok" for r in p):
                paired.append(dict(case=case,seed=seed,native=p[0]["top1_success"],mixed=p[1]["top1_success"],
                                   native5=p[0]["top5_success"],mixed5=p[1]["top5_success"]))
    write_rows(a.output/"pose_pairs.csv",paired)
    print(json.dumps(dict(planned_pairs=len(case_ids)*len(SEEDS),scorable_pairs=len(paired),
          native=sum(p["native"] for p in paired),mixed=sum(p["mixed"] for p in paired),
          native5=sum(p["native5"] for p in paired),mixed5=sum(p["mixed5"] for p in paired)),indent=2))
    if not paired:
        raise RuntimeError("no scorable pose pairs")

def write_rows(path, rows):
    fields = sorted({k for row in rows for k in row})
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def run_affinity(a):
    rows = []
    targets = AFFINITY if a.subset == 'all' else {a.subset: AFFINITY[a.subset]}
    for target, anchor in targets.items():
        data = a.data / "data" / target
        folder = a.output / target
        folder.mkdir(exist_ok=True,parents=True)
        meta = yaml.safe_load((data/"00_data"/"ligands.yml").read_text())
        ref = [m for m in Chem.SDMolSupplier(str(data/"02_ligands"/"ligands.sdf"),removeHs=False)
               if m is not None and m.GetProp("_Name")==anchor]
        if len(ref)!=1:
            raise ValueError("missing fixed pocket anchor")
        center=coords(ref[0]).mean(axis=0)
        pdb=data/"01_protein"/"crd"/"protein.pdb"
        try:
            receptor=prepare_receptor(pdb,folder)
            prep_error=None
        except Exception as e:
            prep_error=repr(e)
        for name, info in sorted(meta.items()):
            m=info["measurement"]
            if m["type"].lower() not in ("ki","kd"):
                continue
            unit={"nM":1e-9,"uM":1e-6,"mM":1e-3,"M":1}.get(m["unit"])
            value=m["value"]
            # Reject qualified/non-positive values; underlying paper qualifiers
            # and assay conditions still need independent adjudication.
            table=re.match(r"table\s+(\w+)",m.get("comment",""),re.I)
            if unit is None or not isinstance(value,(int,float)) or value<=0 or not table:
                continue
            label=dict(target=target,ligand=name,assay=str(m["doi"])+" table "+table.group(1),
                       type=m["type"],ki_M=value*unit,error_M=m["error"]*unit if m.get("error",-1)>=0 else None)
            d=folder/name
            d.mkdir(exist_ok=True)
            ligand=d/"ligand.pdbqt"
            try:
                if prep_error:
                    raise ValueError(prep_error)
                prepare_ligand(info["smiles"],ligand)
                error=None
            except Exception as e:
                error=repr(e)
            for seed, arm in itertools.product(SEEDS, ARMS):
                row=dict(label,seed=seed,arm=arm,status="ok")
                try:
                    if error:
                        raise ValueError(error)
                    _, _, vals=docking(a.binary,receptor,ligand,center,d,seed,arm,a.exhaustiveness,a.cap)
                    row.update(vals)
                except Exception as e:
                    row.update(status="failed",error=repr(e))
                rows.append(row)
            write_rows(a.output/"affinity_rows.csv",rows)
    stats=[]
    for group in sorted({(r["target"],r["assay"],r["type"],r["seed"]) for r in rows}):
        subset=[r for r in rows if (r["target"],r["assay"],r["type"],r["seed"])==group]
        pairs=[]
        for name in sorted({r["ligand"] for r in subset}):
            p=[r for r in subset if r["ligand"]==name and r["status"]=="ok"]
            if len(p)==2:
                pairs.append((p[0],p[1]))
        row=dict(target=group[0],assay=group[1],type=group[2],seed=group[3],
                 n_expected=len({r["ligand"] for r in subset}),n_paired=len(pairs))
        for i,arm in enumerate(ARMS):
            x=[p[i] for p in pairs]
            if len(x)<3: continue
            rho=float(spearmanr([r["ki_M"] for r in x],[r["score"] for r in x]).statistic)
            row[arm+"_rho"]=rho
            n=correct=ties=large_n=large_correct=known_n=known_correct=0
            for u,v in itertools.combinations(x,2):
                if u["ki_M"]==v["ki_M"]: continue
                n+=1
                d=(u["ki_M"]-v["ki_M"])*(u["score"]-v["score"])
                correct+=d>0
                ties+=d==0
                if max(u["ki_M"],v["ki_M"])/min(u["ki_M"],v["ki_M"])>=3:
                    large_n+=1
                    large_correct+=d>0
                if u["error_M"] is not None and v["error_M"] is not None and abs(u["ki_M"]-v["ki_M"])>u["error_M"]+v["error_M"]:
                    known_n+=1
                    known_correct+=d>0
            row.update({arm+"_pairs":n,arm+"_correct":correct,arm+"_ties":ties,
                        arm+"_3fold_pairs":large_n,arm+"_3fold_correct":large_correct,
                        arm+"_nonoverlap_reported_error_pairs":known_n,arm+"_nonoverlap_reported_error_correct":known_correct})
        stats.append(row)
    write_rows(a.output/"affinity_group_stats.csv",stats)
    print("Exploratory source/table-stratified Ki/Kd metrics (not assay-adjudicated).")
    if not any(r["n_paired"]>=3 for r in stats):
        raise RuntimeError("no scorable affinity series")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("binary",type=Path)
    p.add_argument("mode",choices=("pose","affinity"))
    p.add_argument("data",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("--exhaustiveness",type=int,default=4)
    p.add_argument("--cap",type=int,default=3000)
    p.add_argument("--subset",default="all",help="pose shard 0/1/2 or affinity target")
    a=p.parse_args()
    a.binary=a.binary.resolve()
    a.data=a.data.resolve()
    a.output=a.output.resolve()
    a.output.mkdir(parents=True,exist_ok=True)
    (run_pose if a.mode=="pose" else run_affinity)(a)

if __name__=="__main__":
    main()
