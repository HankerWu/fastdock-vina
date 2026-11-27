"""Summarize measured affinity within one assay; verify fixed parent coordinates."""

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import re
import subprocess
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.stats import kendalltau, spearmanr

from assay import MEASURED, OLD, read, run_panel

RT = 0.00198720425864083 * 298.0
ARMS = run_panel.ARMS


def report(args):
    work = args.work
    old = args.original
    out = args.out
    out.mkdir(exist_ok=True, parents=True)
    lookup, measurements = {}, {r["pdb"]: r for r in MEASURED}
    cells, tops = [], {}
    for source in (old, work):
        for row in read(source / "cells.csv"):
            if row["target"] == "FABP4":
                cells.append(row)
                lookup[row["key"]] = source
        for row in read(source / "modes.csv"):
            if row["key"] in lookup and row["mode"] == "1":
                tops[row["key"]] = float(row["score"])
    assert len(tops) == len(cells) == len(MEASURED) * 2 * 2 * 3
    details = []
    for row in cells:
        r = dict(row)
        r["ligand"] = measurements[row["query"]]["id"]
        r["top1_score"] = tops[row["key"]]
        details.append(r)
    run_panel.paired.write(out / "all_cells.csv", details)
    mols = {
        p: Chem.MolFromSmiles(
            json.loads((work / p / "preparation.json").read_text())["canonical_smiles"]
        )
        for p in measurements
    }
    ligands, pairs, ranks, costs = [], [], [], []
    for arm in ARMS:
        allcells = [r for r in details if r["arm"] == arm]
        costs.append(
            dict(
                arm=arm,
                predictions=len(allcells),
                calls=sum(
                    sum(
                        int(r[k])
                        for k in ("hunt_calls", "second_calls", "explicit_calls")
                    )
                    for r in allcells
                ),
                evals=sum(
                    sum(
                        int(r[k])
                        for k in ("hunt_evals", "second_evals", "explicit_evals")
                    )
                    for r in allcells
                ),
                cpu_s=sum(float(r["cpu_s"]) for r in allcells),
                wall_s=sum(float(r["wall_s"]) for r in allcells),
                largest_RSS_KiB=max(int(r["maxrss_KiB"]) for r in allcells),
            )
        )
        for rec in ("both", *OLD):
            group = [r for r in allcells if rec == "both" or r["receiver"] == rec]
            scores = {
                p: np.mean([r["top1_score"] for r in group if r["query"] == p])
                for p in measurements
            }
            for p, l in measurements.items():
                ligands.append(
                    dict(
                        arm=arm,
                        receiver=rec,
                        id=l["id"],
                        pdb=p,
                        Kd_uM=l["Kd_uM"],
                        true_G_kcal=RT * math.log(float(l["Kd_uM"]) * 1e-6),
                        mean_top1_score=scores[p],
                        N_heavy=mols[p].GetNumHeavyAtoms(),
                        charge=Chem.GetFormalCharge(mols[p]),
                    )
                )
            this = []
            for a, b in itertools.combinations(
                sorted(measurements, key=lambda p: float(measurements[p]["Kd_uM"])), 2
            ):
                ratio = float(measurements[b]["Kd_uM"]) / float(
                    measurements[a]["Kd_uM"]
                )
                exp = RT * math.log(ratio)
                pred = scores[b] - scores[a]
                ok = (1 if pred > 0 else 0) if pred != 0 else 0.5
                this.append(
                    dict(
                        arm=arm,
                        receiver=rec,
                        strong=measurements[a]["id"],
                        weak=measurements[b]["id"],
                        Kd_ratio=ratio,
                        measured_delta=exp,
                        pred_delta=pred,
                        equal_N_charge=(
                            mols[a].GetNumHeavyAtoms(),
                            Chem.GetFormalCharge(mols[a]),
                        )
                        == (mols[b].GetNumHeavyAtoms(), Chem.GetFormalCharge(mols[b])),
                        correct=ok,
                    )
                )
            pairs.extend(this)
            for subset, lst in (
                ("all", this),
                ("3fold", [p for p in this if p["Kd_ratio"] >= 3]),
                ("equal_N_charge", [p for p in this if p["equal_N_charge"]]),
            ):
                ranks.append(
                    dict(
                        arm=arm,
                        receiver=rec,
                        subset=subset,
                        pairs=len(lst),
                        correct=sum(p["correct"] for p in lst),
                        mean_abs_delta_error_kcal=np.mean(
                            [abs(p["pred_delta"] - p["measured_delta"]) for p in lst]
                        ),
                        spearman=(
                            spearmanr(
                                [float(l["Kd_uM"]) for l in measurements.values()],
                                list(scores.values()),
                            ).statistic
                            if subset == "all"
                            else ""
                        ),
                        kendall_b=(
                            kendalltau(
                                [float(l["Kd_uM"]) for l in measurements.values()],
                                list(scores.values()),
                            ).statistic
                            if subset == "all"
                            else ""
                        ),
                    )
                )
    for file, rows in (
        ("ligands", ligands),
        ("pairs", pairs),
        ("ranks", ranks),
        ("costs", costs),
    ):
        run_panel.paired.write(out / (file + ".csv"), rows)
    if args.fixed:
        score_fixed(args, details, lookup)


def score_fixed(args, details, lookup):
    d = args.out / "fixed"
    d.mkdir(exist_ok=True)
    os.sched_setaffinity(0, {args.cpu})
    rows = []
    for row in details:
        if row["arm"] != "parent":
            continue
        key = row["key"]
        raw = (lookup[key] / "runs" / (key + ".pdbqt")).read_text().split("ENDMDL")[0]
        ligand = d / (key + ".pdbqt")
        ligand.write_text(
            "\n".join(l for l in raw.splitlines() if l and not l.startswith("MODEL"))
            + "\n"
        )
        unbound = float(re.search(r"REMARK UNBOUND:\s*(-?[\d.]+)", raw).group(1))
        common = [
            "--score_only",
            "--cpu",
            "1",
            "--seed",
            "1",
            "--receptor",
            str(args.work / row["receiver"] / "receptor.pdbqt"),
            "--ligand",
            str(ligand),
            "--unbound_energy=" + str(unbound),
        ]
        for axis, c in zip(
            "xyz",
            json.loads((args.work / row["receiver"] / "preparation.json").read_text())[
                "center"
            ],
        ):
            common.extend(["--center_" + axis, str(c), "--size_" + axis, "22.5"])
        for arm in ARMS:
            proc = subprocess.run(
                [str(getattr(args, arm).resolve()), *common],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=True,
            )
            (d / (key + "_" + arm + ".log")).write_text(proc.stdout)
            nums = [
                float(
                    re.search(
                        re.escape(label) + r"\s*:\s*([-\d.]+)", proc.stdout
                    ).group(1)
                )
                for label in (
                    "Estimated Free Energy of Binding",
                    "Final Intermolecular Energy",
                    "Final Total Internal Energy",
                    "Torsional Free Energy",
                    "Unbound System's Energy",
                )
            ]
            rows.append(
                dict(
                    key=key,
                    arm=arm,
                    score=nums[0],
                    inter=nums[1],
                    intra=nums[2],
                    torsional=nums[3],
                    unbound=nums[4],
                    pose_sha256=hashlib.sha256(ligand.read_bytes()).hexdigest(),
                )
            )
    run_panel.paired.write(args.out / "fixed_top1.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("work", "original", "out"):
        p.add_argument("--" + key, type=Path, required=True)
    p.add_argument("--fixed", action="store_true")
    p.add_argument("--cpu", type=int, default=4)
    for arm in ARMS:
        p.add_argument("--" + arm, type=Path)
    report(p.parse_args())
