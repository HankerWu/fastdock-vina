"""Generate fixed validation pools, then evaluate registered rankings on those exact poses."""
import argparse
import csv
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH / "frozen_core/hydroxyl_partition"))
from scorer import Receptor, score_pose

BUDGETS = (3000, 12000)
SEEDS = (43103, 43117)


def cases(work):
    systems = json.loads((work / "sources.json").read_text())["systems"]
    for pdb, record in systems.items():
        if record["receptor_status"] == "ready":
            yield pdb, json.loads((work / pdb / "preparation.json").read_text())


def dump(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def generate(work, pro, upstream):
    where = work / "runs"
    where.mkdir(exist_ok=True)
    costs = []
    for pdb, meta in cases(work):
        receptor, ligand = (work / pdb / (n + ".pdbqt") for n in ("receptor", "ligand"))
        for cap, seed, (name, binary) in product(BUDGETS, SEEDS, (("pro", pro), ("upstream", upstream))):
            stem = where / f"{pdb}_{cap}_{seed}_{name}"
            cli = [str(binary.resolve()), "--receptor", str(receptor), "--ligand", str(ligand),
                   "--cpu", "1", "--exhaustiveness", "4", "--max_evals", str(cap),
                   "--num_modes", "9", "--energy_range", "100", "--seed", str(seed),
                   "--out", str(stem.with_suffix(".pdbqt"))]
            if name == "pro":
                cli.append("--search_counters")
            for axis, center in zip("xyz", meta["center"]):
                cli.extend(["--size_" + axis, "24", "--center_" + axis, str(center)])
            log = stem.with_suffix(".log")
            cost = stem.with_suffix(".json")
            subprocess.run([sys.executable, "-S", str(BENCH / "rigid_flip/process_cost.py"), log, cost, *cli], check=True)
            totals = re.search(r"Search counters: (.+)", log.read_text())
            costs.append(dict(pdb=pdb, target=meta["target"], cap=cap, seed=seed, engine=name,
                              search_counters=totals.group(1).strip() if totals else "",
                              **json.loads(cost.read_text())))
    dump(work / "docking_costs.csv", costs)
    binaries = {label: {"version": subprocess.check_output([path, "--version"], text=True).strip(),
                        "sha256": sha256(path.read_bytes()).hexdigest()}
                for label, path in (("pro", pro), ("upstream", upstream))}
    (work / "binaries.json").write_text(json.dumps(binaries, indent=2) + "\n")


def evaluate(work, rotor_probe, with_bonded=False):
    if with_bonded:
        sys.path.insert(0, str(BENCH / "frozen_core/hydroxyl_evidence"))
        from bonded import add_bonded_prior
    all_modes, rank_cells, costs = [], [], []
    for pdb, meta in cases(work):
        ligand = work / pdb / "ligand.pdbqt"
        protein_text = (work / pdb / "receptor.pdbqt").read_text()
        ligand_text = ligand.read_text()
        prior_start = time.perf_counter()
        receptor = Receptor(protein_text)
        prior_s = time.perf_counter() - prior_start
        bonded_start = time.perf_counter()
        bonded = add_bonded_prior(receptor) if with_bonded else None
        bonded_s = time.perf_counter() - bonded_start if with_bonded else 0.
        readout = subprocess.check_output([rotor_probe, ligand], text=True)
        torsions = float(next(csv.DictReader(readout.splitlines()))["num_tors"])
        divisor = 1. + .05846 * torsions
        for cap, seed, engine in product(BUDGETS, SEEDS, ("pro", "upstream")):
            key = dict(pdb=pdb, target=meta["target"], cap=cap, seed=seed, engine=engine)
            path = work / "runs" / f"{pdb}_{cap}_{seed}_{engine}.pdbqt"
            modes = []
            scoring_s = 0.
            bonded_scoring_s = 0.
            for i, block in enumerate(path.read_text().split("ENDMDL")[:-1], 1):
                native = float(next(s for s in block.splitlines() if s.startswith("REMARK VINA RESULT")).split()[3])
                xyz = {int(s[6:11]): [float(s[q:q+8]) for q in (30, 38, 46)]
                       for s in block.splitlines() if s.startswith(("ATOM  ", "HETATM"))}
                heavy = np.array([xyz[s] for s in meta["serial"]])
                rmsd = np.sqrt(((heavy - np.array(meta["mapped"]))**2).sum(axis=2).mean(axis=1).min())
                mode = dict(**key, native_rank=i, native=native, rmsd=float(rmsd))
                if engine == "pro":
                    start = time.perf_counter()
                    parts = score_pose(receptor, ligand_text, block)
                    scoring_s += time.perf_counter() - start
                    mode.update(parts)
                    mode["num_tors"] = torsions
                    for label in ("fixed", "minimum", "uniform", "ensemble", "rotor_only", "nonrotor"):
                        mode[label] = native + parts["raw_" + label] / divisor
                    if bonded is not None:
                        start = time.perf_counter()
                        mode["raw_bonded"] = score_pose(bonded, ligand_text, block)["raw_ensemble"]
                        mode["bonded"] = native + mode["raw_bonded"] / divisor
                        bonded_scoring_s += time.perf_counter() - start
                modes.append(mode)
            arms = ["native"]
            if engine == "pro":
                arms.extend(["fixed", "minimum", "uniform", "ensemble", "rotor_only", "nonrotor"])
                if bonded is not None:
                    arms.append("bonded")
                costs.append(dict(**key, modes=len(modes), receptor_once_s=prior_s, score_total_s=scoring_s,
                                  bonded_additional_receptor_s=bonded_s, bonded_score_total_s=bonded_scoring_s))
            for arm in arms:
                ordered = modes if arm == "native" else sorted(modes, key=lambda m: (m[arm], m["native_rank"]))
                for rank, m in enumerate(ordered, 1):
                    m["rank_" + arm] = rank
                rms = [m["rmsd"] for m in ordered]
                rank_cells.append(dict(**key, arm=arm, modes=len(rms), top1_mode=ordered[0]["native_rank"],
                                       rmsd_top1=rms[0], rmsd_best5=min(rms[:5]), rmsd_oracle=min(rms),
                                       success1=int(rms[0]<=2), success5=int(min(rms[:5])<=2), oracle=int(min(rms)<=2),
                                       top1_gap=(ordered[1][arm] - ordered[0][arm]) if len(rms)>1 else ""))
            all_modes.extend(modes)
    # Keep a rectangular table, including upstream rows with native rank only.
    fields = dict.fromkeys(k for m in all_modes for k in m)
    dump(work / "modes.csv", [{k:m.get(k, "") for k in fields} for m in all_modes])
    dump(work / "cells.csv", rank_cells)
    dump(work / "scoring_costs.csv", costs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--pro", type=Path)
    p.add_argument("--upstream", type=Path)
    p.add_argument("--rotor-probe", type=Path)
    p.add_argument("--bonded", action="store_true", help="Include the PPAR challenge bonded-prior scoring arm")
    args = p.parse_args()
    if bool(args.pro) != bool(args.upstream):
        p.error("Supply both native engines to generate their pools")
    if args.pro:
        generate(args.work.resolve(), args.pro, args.upstream)
    if args.rotor_probe:
        evaluate(args.work.resolve(), args.rotor_probe.resolve(), args.bonded)
