"""Run frozen paired predictions. Each search owns one core; seed never selects a pose."""

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from time import monotonic

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "cubic_grid/flex_replay"))
from run_flex import reference
from prepare import MEASUREMENTS

SEEDS = (62483, 87911, 120997)
HISTORICAL = {"FABP4": ("2NNQ", "3FR2"), "RBP4": ("1RBP", "3FMZ")}


def csvout(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def cases(prepared, solutes=False):
    load = lambda p: json.loads((prepared / p / "metadata.json").read_text())
    for row in json.loads((prepared / "manifest.json").read_text()):
        if row["status"] == "ready":
            p = row["pdb"]
            m = load(p)
            yield (
                ("crystal_additive" if solutes else "new"),
                row["target"],
                p,
                p,
                prepared / p / "ligand.pdbqt",
                m["serial"],
                m["mapped"],
            )
    if solutes:
        return
    for target, group in HISTORICAL.items():
        for q, r in itertools.product(group, repeat=2):
            qm, rm = load(q), load(r)
            coords = qm["mapped"] if q == r else reference(qm, rm)
            yield (
                "historical",
                target,
                q,
                r,
                prepared / q / "ligand.pdbqt",
                qm["serial"],
                coords,
            )
    for label in csv.DictReader(MEASUREMENTS.open()):
        for r in HISTORICAL["FABP4"]:
            q = label["id"]
            yield (
                "assay",
                "FABP4",
                q,
                r,
                prepared / "assay" / (q + ".pdbqt"),
                None,
                None,
            )


def measure(path, serial, reference_xyz):
    out = []
    for block in path.read_text().split("ENDMDL"):
        ligand = block.split("BEGIN_RES")[0].splitlines()
        scores = [
            float(line.split()[3])
            for line in ligand
            if line.startswith("REMARK VINA RESULT:")
        ]
        if not scores:
            continue
        atoms = {
            int(l[6:11]): [float(l[i : i + 8]) for i in (30, 38, 46)]
            for l in ligand
            if l.startswith(("ATOM  ", "HETATM"))
        }
        row = dict(rank=len(out) + 1, score=scores[0])
        if serial is not None:
            observed = np.asarray([atoms[i] for i in serial])
            distance = np.asarray(reference_xyz) - observed
            row["rmsd"] = float(
                np.sqrt(np.mean(np.sum(distance**2, axis=2), axis=1)).min()
            )
        out.append(row)
    return out


def execute(job):
    prepared, raw, arm, binary, case, seed, budget = job
    cohort, target, q, r, lig, serial, expected = case
    meta = json.loads((prepared / r / "metadata.json").read_text())
    key = dict(cohort=cohort, target=target, query=q, receiver=r, seed=seed, arm=arm)
    stem = raw / f"{cohort}-{q}-{r}-{seed}-{arm}"
    pose, log = stem.with_suffix(".pdbqt"), stem.with_suffix(".log")
    cmd = [
        binary,
        "--receptor",
        prepared / r / "partition_rigid.pdbqt",
        "--flex",
        prepared / r / "partition_flex.pdbqt",
        "--ligand",
        lig,
        "--cpu",
        "1",
        "--exhaustiveness",
        "1",
        "--max_evals",
        str(budget),
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--seed",
        str(seed),
        "--out",
        pose,
    ]
    for axis, value in zip("xyz", meta["center"]):
        cmd += ["--center_" + axis, str(value), "--size_" + axis, "24"]
    if arm != "official":
        cmd.append("--search_counters")
    start = monotonic()
    with log.open("w") as f:
        p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
        _, st, cost = os.wait4(p.pid, 0)
        p.returncode = os.waitstatus_to_exitcode(st)
    key.update(
        wall_s=monotonic() - start,
        cpu_s=cost.ru_utime + cost.ru_stime,
        peak_kib=cost.ru_maxrss,
        exit=p.returncode,
    )
    counts = re.findall(
        r"(hunt_evals|second_evals|explicit_evals)=(\d+)", log.read_text()
    )
    key.update({k: int(v) for k, v in counts})
    modes = (
        measure(pose, serial, expected) if not p.returncode and pose.exists() else []
    )
    key.update(reported=len(modes), top1_score=modes[0]["score"] if modes else None)
    if modes and serial is not None:
        rms = [r["rmsd"] for r in modes]
        key.update(
            top1=rms[0],
            best5=min(rms[:5]),
            best_reported=min(rms),
            regret=rms[0] - min(rms),
        )
    keys = {k: key[k] for k in ("cohort", "target", "query", "receiver", "seed", "arm")}
    return key, [dict(**keys, **x) for x in modes]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("prepared", "out", "baseline", "fixed", "official"):
        p.add_argument(k, type=Path)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--solutes", action="store_true")
    p.add_argument(
        "--extended-budget",
        action="store_true",
        help="24000 budget paired Pro structure runs",
    )
    a = p.parse_args()
    w, o = a.prepared.resolve(), a.out.resolve()
    raw = o / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    arms = {k: getattr(a, k).resolve() for k in ("baseline", "fixed", "official")}
    (o / "binaries.json").write_text(
        json.dumps(
            {k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in arms.items()},
            indent=2,
        )
        + "\n"
    )
    jobs = []
    for i, (case, seed) in enumerate(itertools.product(cases(w, a.solutes), SEEDS)):
        if a.extended_budget and case[0] == "assay":
            continue
        order = (
            list(arms)
            if case[0] != "assay" and not a.extended_budget
            else ["baseline", "fixed"]
        )
        order = order[i % len(order) :] + order[: i % len(order)]
        jobs.extend(
            (w, raw, k, arms[k], case, seed, 24000 if a.extended_budget else 12000)
            for k in order
        )
    with ThreadPoolExecutor(a.jobs) as workers:
        result = list(workers.map(execute, jobs))
    csvout(o / "cases.csv", [x[0] for x in result])
    csvout(o / "modes.csv", [m for x in result for m in x[1]])


if __name__ == "__main__":
    main()
