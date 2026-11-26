"""Run paired hCA II full searches or rescore their output, including failed jobs."""

import argparse
import concurrent.futures
import csv
import itertools
import json
import re
import subprocess
import time
from pathlib import Path

import numpy as np
from prep import CASES, RECS
from rdkit import Chem
from rdkit.Chem import rdMolAlign

S = [42, 1729]
ARMS = ["native", "smooth"]


def tasks(cognate=False):
    pairs = (
        ((p, p) for p in CASES if p not in RECS)
        if cognate
        else itertools.product(RECS, CASES)
    )
    return [
        (r, l, s, a)
        for (r, l), (s, a) in itertools.product(pairs, itertools.product(S, ARMS))
    ]


def job(task, binary, prepared, out):
    rec, lig, seed, arm = task
    dest = out / f"{rec}_{lig}_{seed}_{arm}"
    qt = dest.with_suffix(".pdbqt")
    log = dest.with_suffix(".log")
    cmd = [
        str(binary.resolve()),
        "--receptor",
        str(prepared / f"{rec}_receptor.pdbqt"),
        "--ligand",
        str(prepared / f"{lig}_ligand.pdbqt"),
        "--config",
        str(prepared / f"{rec}_box.txt"),
        "--seed",
        str(seed),
        "--cpu",
        "1",
        "--scoring",
        "vina",
        "--exhaustiveness",
        "4",
        "--max_evals",
        "50000",
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--search_counters",
        "--out",
        str(qt),
    ]
    if arm == "smooth":
        cmd += ["--smooth_cutoff"]
    ts = time.monotonic()
    with log.open("w") as f:
        c = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if c.returncode:
        raise ValueError(f"Failed: {log}: " + log.read_text()[-900:])
    print(rec, lig, seed, arm, round(time.monotonic() - ts, 1), flush=True)


def evaluate(prepared, out, cognate=False):
    data = json.loads((prepared / "metadata.json").read_text())
    rows = []
    mods = []
    for rec, lig, seed, arm in tasks(cognate):
        dest = out / f"{rec}_{lig}_{seed}_{arm}"
        newrows = []
        conf = data[lig]["references"][rec]["coords"]
        d = data[lig]
        ref = Chem.SDMolSupplier(str(prepared / f"{lig}_graph.sdf"), removeHs=True)[0]
        rc = ref.GetConformer()
        for i, pt in enumerate(conf):
            rc.SetAtomPosition(i, pt)
        pr = Chem.Mol(ref)
        pc = pr.GetConformer()
        for chunk in dest.with_suffix(".pdbqt").read_text().split("ENDMDL"):
            en = re.search(r"^REMARK VINA RESULT:\s+(\S+)", chunk, re.M)
            if not en:
                continue
            xyz = {
                int(l[6:11]): [float(l[30 + i : 38 + i]) for i in (0, 8, 16)]
                for l in chunk.splitlines()
                if l.startswith(("ATOM", "HETATM"))
            }
            for i, ser in enumerate(d["qt_serial"]):
                pc.SetAtomPosition(i, xyz[ser])
            rms = rdMolAlign.CalcRMS(
                pr, ref, maxMatches=1000000, symmetrizeConjugatedTerminalGroups=True
            )
            newrows.append(
                dict(
                    rec=rec,
                    lig=lig,
                    Ki_nM=d["Ki_nM"],
                    seed=seed,
                    arm=arm,
                    model=len(newrows) + 1,
                    score=float(en[1]),
                    rmsd=round(rms, 6),
                )
            )
        assert newrows
        mods.extend(newrows)
        log = dest.with_suffix(".log").read_text()
        ev = {
            key: sum(map(int, re.findall(r"\b" + key + r"=(\d+)", log)))
            for key in ["hunt_evals", "second_evals", "explicit_evals"]
        }
        rows.append(
            dict(
                **newrows[0],
                best_rmsd=min(z["rmsd"] for z in newrows),
                modes=len(newrows),
                **ev,
            )
        )
    for filename, l in [("all_modes.csv", mods), ("rows.csv", rows)]:
        with (out / filename).open("w") as f:
            w = csv.DictWriter(f, fieldnames=list(l[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(l)
    if cognate:
        return
    pairs = []
    groups = []
    for rec, seed, arm in itertools.product(RECS, S, ARMS):
        rr = [r for r in rows if (r["rec"], r["seed"], r["arm"]) == (rec, seed, arm)]
        correct = 0
        for l, r in itertools.combinations(sorted(rr, key=lambda r: r["Ki_nM"]), 2):
            if r["Ki_nM"] / l["Ki_nM"] < 3:
                continue
            diff = r["score"] - l["score"]
            val = 1 if diff > 0 else 0.5 if diff == 0 else 0
            correct += val
            pairs.append(
                dict(
                    rec=rec,
                    seed=seed,
                    arm=arm,
                    strong=l["lig"],
                    weak=r["lig"],
                    Ki_ratio=r["Ki_nM"] / l["Ki_nM"],
                    score_weak_minus_strong=round(diff, 4),
                    credit=val,
                )
            )
        groups.append(
            dict(
                rec=rec,
                seed=seed,
                arm=arm,
                pose_top1=sum(r["rmsd"] <= 2 for r in rr),
                pose_oracle=sum(r["best_rmsd"] <= 2 for r in rr),
                n=len(rr),
                median_top1_rmsd=float(np.median([r["rmsd"] for r in rr])),
                pairs_correct=correct,
                pairs=sum(
                    (p["rec"], p["seed"], p["arm"]) == (rec, seed, arm) for p in pairs
                ),
            )
        )
    for filename, l in [("pairs.csv", pairs), ("groups.csv", groups)]:
        with (out / filename).open("w") as f:
            w = csv.DictWriter(f, fieldnames=list(l[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(l)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--evaluate-only", action="store_true")
    parser.add_argument("--extra-cognate", action="store_true")
    args = parser.parse_args()
    prepared = args.work / "prepared"
    out = args.output or args.work / ("cognate" if args.extra_cognate else "runs")
    if not args.evaluate_only:
        if args.binary is None:
            parser.error("--binary is required for docking")
        out.mkdir(parents=True, exist_ok=False)
        with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
            futures = [
                pool.submit(job, t, args.binary, prepared, out)
                for t in tasks(args.extra_cognate)
            ]
            for future in concurrent.futures.as_completed(futures):
                future.result()
    evaluate(prepared, out, args.extra_cognate)


if __name__ == "__main__":
    main()
