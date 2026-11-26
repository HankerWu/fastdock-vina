#!/usr/bin/env python3
"""Prospective small cognate panel for geometric BFGS initialization.

Fixed public PDB/CCD cases; no crystal-derived starting conformations.
Reject dry single-chain preparations missing a nonwater atom within 8 A.
The panel tests generalization across targets, not cross-docking or affinity.
"""
import argparse
import csv
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmarks/paired_crossdocking"))
import run as base

CASES = (
    ("5SAK", "ZRY", "endothiapepsin"),
    ("5SB2", "1K2", "DDR1"),
    ("5SD5", "HWI", "DHFR_5SD5"),
    ("5SIS", "JSM", "PDE10A"),
    ("6M2B", "EZO", "6M2B"),
    ("6M73", "FNR", "lactate_oxidase"),
    ("6T88", "MWQ", "urocanate_reductase"),
    ("6XBO", "5MC", "CMP_sialic_acid_transporter"),
    ("6XCT", "478", "pepsin_A"),
    ("6XG5", "TOP", "E_coli_DHFR"),
)
SEEDS = (63419, 63421)
ARMS = ("armijo", "geometry")


def prepare(work, pdb, ligand, target):
    out = work / pdb
    out.mkdir(parents=True, exist_ok=True)
    source, ccd = out / (pdb + ".pdb"), out / (ligand + ".cif")
    row = dict(pdb=pdb, ligand=ligand, target=target)
    row["pdb_sha256"] = base.fetch(source, f"https://files.rcsb.org/download/{pdb}.pdb")
    row["ccd_sha256"] = base.fetch(ccd, f"https://files.rcsb.org/ligands/download/{ligand}.cif")
    lines = source.read_text().splitlines()
    ligands = [l for l in lines if l.startswith("HETATM") and l[17:20].strip() == ligand
               and l[76:78].strip() not in ("H", "D")]
    if not ligands:
        raise ValueError("no deposited ligand")
    instances = sorted({l[21:27] for l in ligands})
    # Deterministically select the lexicographically first deposited instance.
    instance = instances[0]
    instance_atoms = [l for l in ligands if l[21:27] == instance]
    # Named conformers may start at B, C or D. Match dep_records and
    # include blank atoms shared by every named conformer.
    labels = sorted({l[16] for l in instance_atoms} - {" "})
    first_alt = labels[0] if labels else " "
    obs = np.array([base.coords(l) for l in instance_atoms if l[16] in (" ", first_alt)])
    protein = [l for l in lines if l.startswith("ATOM  ") and l[76:78].strip() not in ("H", "D")]
    chains = sorted({l[21] for l in protein})
    distance = lambda l: float(np.linalg.norm(obs - np.asarray(base.coords(l)), axis=1).min())
    chain = min(chains, key=lambda c: min(distance(l) for l in protein if l[21] == c))
    omit = [l for l in lines if l.startswith(("ATOM  ", "HETATM"))
            and l[76:78].strip() not in ("H", "D")
            and not (l.startswith("ATOM  ") and l[21] == chain)
            and not (l.startswith("HETATM") and l[17:20].strip() == ligand and l[21:27] == instance)
            and l[17:20].strip() not in ("HOH", "DOD", "WAT")]
    nearest = min((distance(l) for l in omit), default=math.inf)
    row.update(chain=chain, ligand_instance=instance,
               nearest_omitted_nonwater_A=None if not math.isfinite(nearest) else nearest)
    (out / "source.json").write_text(json.dumps(row, indent=2))
    if nearest < 8:
        raise ValueError(f"omitted nonwater atom within 8 A: {nearest:.3f}")
    # Reuse the established chemistry mapping and preparation, remapping the
    # one selected receptor chain and one ligand instance to chain A.
    selected = []
    for l in lines:
        if l.startswith("ATOM  ") and l[21] == chain:
            selected.append(l[:21] + "A" + l[22:])
        elif l.startswith("HETATM") and l[17:20].strip() == ligand and l[21:27] == instance:
            selected.append(l[:21] + "A" + l[22:])
    refs = base.dep_records(selected, ligand)
    if not refs[0]:
        raise ValueError("no complete ligand reference")
    row.update(base.ligand(ccd, refs, out / "ligand.pdbqt"))
    row.update(base.receptor(selected, out, pdb, obs))
    repaired = row.get("closest_repaired_ligand_distance")
    if repaired is not None and repaired < 6:
        raise ValueError(f"repaired heavy atom within 6 A: {repaired:.3f}")
    row["center"] = obs.mean(axis=0).tolist()
    row["status"] = "ready"
    (out / "metadata.json").write_text(json.dumps(row, indent=2))
    return row


def dock(work, binary, meta, arm, seed):
    out = work / meta["pdb"]
    path = out / f"{seed}-{arm}.pdbqt"
    cmd = [str(binary), "--receptor", str(out/"receptor.pdbqt"), "--ligand", str(out/"ligand.pdbqt"),
           "--pose_selection", "late", "--refinement_optimizer", arm, "--exhaustiveness", "3",
           "--max_evals", "12000", "--num_modes", "9", "--energy_range", "100",
           "--search_counters", "--cpu", "1", "--seed", str(seed), "--out", str(path)]
    for ax, center in zip("xyz", meta["center"]):
        cmd += ["--center_"+ax, str(center), "--size_"+ax, "24"]
    start = time.monotonic()
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=300)
    wall = time.monotonic()-start
    path.with_suffix(".log").write_text(res.stdout)
    res.check_returncode()
    reference = np.asarray(meta["mapped"])
    values = []
    for block in path.read_text().split("ENDMDL"):
        atoms = {int(l[6:11]): base.coords(l) for l in block.splitlines()
                 if l.startswith(("ATOM  ", "HETATM"))}
        if not atoms:
            continue
        xyz = np.asarray([atoms[i] for i in meta["serial"]])
        rms = float(np.sqrt(np.mean(np.sum((reference-xyz)**2, axis=2),axis=1)).min())
        # Independent arithmetic check of mapping-based receptor-frame RMSD.
        scalar = min(math.sqrt(sum(sum((float(x)-float(y))**2 for x,y in zip(a,b))
                                  for a,b in zip(mapping,xyz))/len(xyz))
                     for mapping in reference)
        if abs(rms-scalar) > 1e-9:
            raise ValueError("RMSD recomputation mismatch")
        score = float(next(l.split()[3] for l in block.splitlines() if l.startswith("REMARK VINA RESULT:")))
        values.append(dict(pdb=meta["pdb"],target=meta["target"],arm=arm,seed=seed,
                           rank=len(values)+1,score=score,rmsd=rms,wall_s=wall))
    if not values:
        raise ValueError("no emitted models")
    return values


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("binary",type=Path)
    a = p.parse_args()
    work = a.work.resolve()
    work.mkdir(parents=True,exist_ok=True)
    manifest, modes, cells = [], [], []
    for pdb, ligand, target in CASES:
        try:
            meta = prepare(work,pdb,ligand,target)
        except Exception as exc:
            manifest.append(dict(pdb=pdb,ligand=ligand,target=target,status="excluded",reason=repr(exc)))
            continue
        manifest.append(dict(pdb=pdb,ligand=ligand,target=target,status="ready"))
        for i,seed in enumerate(SEEDS):
            order = ARMS if i%2==0 else ARMS[::-1]
            for arm in order:
                try:
                    rows = dock(work,a.binary.resolve(),meta,arm,seed)
                    modes.extend(rows)
                    cells.append(dict(pdb=pdb,seed=seed,arm=arm,status="ok"))
                except Exception as exc:
                    cells.append(dict(pdb=pdb,seed=seed,arm=arm,status="error",reason=repr(exc)))
    (work/"manifest.json").write_text(json.dumps(manifest,indent=2))
    (work/"cells.json").write_text(json.dumps(cells,indent=2))
    if modes:
        with (work/"modes.csv").open("w",newline="") as stream:
            writer = csv.DictWriter(stream,fieldnames=list(modes[0]))
            writer.writeheader()
            writer.writerows(modes)
    paired = []
    for pdb,_,_ in CASES:
        for seed in SEEDS:
            r = {arm: [row for row in modes if row["pdb"]==pdb and row["seed"]==seed and row["arm"]==arm]
                 for arm in ARMS}
            if all(r.values()):
                paired.append(dict(pdb=pdb,seed=seed,**{arm: dict(top1=r[arm][0]["rmsd"],
                    oracle=min(m["rmsd"] for m in r[arm]),score=r[arm][0]["score"],
                    wall_s=r[arm][0]["wall_s"]) for arm in ARMS}))
    summary = dict(manifest=manifest,paired_cells=len(paired),pairs=paired)
    for arm in ARMS:
        summary[arm] = dict(top1_hits=sum(r[arm]["top1"]<=2 for r in paired),
                            oracle_hits=sum(r[arm]["oracle"]<=2 for r in paired),
                            total_wall_s=sum(r[arm]["wall_s"] for r in paired))
    (work/"summary.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    if len({r["pdb"] for r in paired}) < 3:
        raise SystemExit("Insufficient eligible structures; no panel conclusion")


if __name__=="__main__":
    main()
