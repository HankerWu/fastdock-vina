"""Run the predeclared 36 blind CLI jobs; analyze their actual released pools."""

import argparse
import csv
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "mutation_scan"))
import run as scan
import quality_cost

spec = importlib.util.spec_from_file_location("independent_metrics", HERE.parent / "mutation_controls/run.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

ARMS = {
    "official": ("official", 25000), "official_double": ("official", 50000),
    "pro": ("base", 25000), "pro_double": ("base", 50000),
    "sham": ("sham", 25000), "pivot": ("pivot", 25000),
}
SEEDS = (81371, 81373)


def totals(log):
    text = log.read_text()
    return {
        key: sum(map(int, re.findall(r"\b" + key + r"=(\d+)", text)))
        for key in ("hunt_calls", "hunt_evals", "second_calls", "second_evals",
                    "explicit_calls", "explicit_evals", "late_prescore")
    }


def dump(path, data):
    path.write_text(json.dumps(data, indent=2) + "\n")


def calculate(args):
    out = args.out
    raw = out / "raw"
    raw.mkdir(parents=True, exist_ok=False)
    counters = out / "official_counter_replay"
    counters.mkdir()
    binaries = {s: args.build / s / "vina" for s in ("base", "sham", "pivot")}
    binaries["official"] = args.official
    if scan.replay.digest(args.official) != audit.OFFICIAL:
        raise ValueError("Expected the official Vina v1.2.7 Linux release")
    systems = scan.config(args.aurora, args.brd9, out)
    frozen = json.loads((HERE.parent / "mutation_scan/results/inputs.json").read_text())
    for key, system in systems.items():
        identity = {k: scan.replay.digest(system[k]) for k in ("ligand", "bonds", "receptor")}
        if any(identity[k] != frozen[key][k] for k in identity):
            raise ValueError(f"Changed input {key}")
    dump(out / "inputs.json", {
        "inputs": {k: frozen[k] for k in systems}, "reference_rmsd": "Receptor frame, graph symmetry",
        "binary_sha256": {k: scan.replay.digest(v) for k, v in dict(binaries, official_counter=args.counter).items()},
        "seeds": SEEDS, "arms": ARMS, "threads": 1, "trajectories": 2, "modes": 9, "energy_range": 100,
    })
    jobs, measured, replay_cost = [], [], []
    for di, (drug, system) in enumerate(systems.items()):
        for si, seed in enumerate(SEEDS):
            order = list(ARMS)
            offset = 2 * di + si
            for arm in order[offset:] + order[:offset]:
                executable, budget = ARMS[arm]
                official = executable == "official"
                stem = f"{drug}_{seed}_{arm}"
                pool = raw / (stem + ".pdbqt")
                log = pool.with_suffix(".log")
                cmd = scan.command(binaries[executable], system, seed, pool, official=official, cap=budget)
                cost = scan.replay.invoke(cmd, log)
                key = dict(drug=drug, seed=seed, arm=arm)
                jobs.append(dict(**key, nominal_evals_per_trajectory=budget, **cost, pool_sha256=scan.replay.digest(pool)))
                if official:
                    other = counters / pool.name
                    counter_cmd = scan.command(args.counter, system, seed, other, official=True, cap=budget)
                    check = scan.replay.invoke(counter_cmd, other.with_suffix(".log"))
                    identical = other.read_bytes() == pool.read_bytes()
                    replay_cost.append(dict(**key, **check, pool_sha256=scan.replay.digest(other),
                                            official_byte_equal=int(identical), **totals(other.with_suffix(".log"))))
                    measurements = totals(other.with_suffix(".log")) if identical else {k: None for k in totals(log)}
                else:
                    measurements = totals(log)
                if measurements["hunt_calls"] == 0:
                    raise ValueError(f"Missing search counters for {key}")
                measured.append(dict(**key, **measurements))
    for name, rows in (("jobs.csv", jobs), ("counts.csv", measured), ("counter_replay.csv", replay_cost)):
        scan.replay.geom.write_csv(out, name, rows)


def analyze(args):
    out = args.out
    systems = scan.config(args.aurora, args.brd9, out)
    independent = audit.inputs(args.aurora, args.brd9, out)
    jobs = list(csv.DictReader((out / "jobs.csv").open()))
    counts = list(csv.DictReader((out / "counts.csv").open()))
    modes, cells = [], []
    worst = 0
    for job, c in zip(jobs, counts, strict=True):
        drug, seed, arm = (job[k] for k in ("drug", "seed", "arm"))
        key = dict(drug=drug, seed=seed, arm=arm)
        assert all(c[k] == job[k] for k in key)
        pool = out / "raw" / f"{drug}_{seed}_{arm}.pdbqt"
        if scan.replay.digest(pool) != job["pool_sha256"]:
            raise ValueError(f"Changed pool {pool}")
        rms = scan.metrics(systems[drug], pool)
        _, ligand, _, _, _, mol, serial, ref = independent[drug]
        checked = audit.measure(pool, ligand, mol, serial, ref)
        for rank, ((r, score), (r_ind, score_ind, drift)) in enumerate(zip(rms, checked, strict=True), 1):
            delta = abs(r - r_ind)
            worst = max(worst, delta)
            assert delta < 1e-4 and score == score_ind
            modes.append(dict(**key, rank=rank, rmsd_A=r, rdkit_independent_A=r_ind,
                              energy=score, max_heavy_bond_drift_A=drift))
        errors = [r for r, _ in rms]
        evals = [c[k] for k in ("hunt_evals", "second_evals", "explicit_evals")]
        cells.append(dict(**key, modes=len(rms), top1=errors[0] if rms else None,
                          best5=min(errors[:5], default=None), oracle=min(errors, default=None),
                          regret=errors[0] - min(errors) if rms else None,
                          evals=sum(map(int, evals)) if all(evals) else None,
                          wall_s=float(job["wall_s"]), cpu_s=float(job["cpu_s"])))
    for name, rows in (("cells.csv", cells), ("modes.csv", modes)):
        scan.replay.geom.write_csv(out, name, rows)
    # Validity, graph-preserving bond lengths, donor/acceptor-independent
    # steric distances, and orientation of every defined heavy stereocentre.
    quality_cost.summary(args.aurora, args.brd9, out)
    distribution = []
    for arm in ARMS:
        selected = [c for c in cells if c["arm"] == arm]
        ok = [c for c in selected if c["modes"]]
        top = np.array([c["top1"] for c in ok])
        distribution.append(dict(arm=arm, cells=len(selected),
                                 **{f"mean_{k}": np.mean([c[k] for c in ok])
                                    for k in ("top1", "best5", "oracle", "regret")},
                                 **dict(zip(("q25", "median", "q75", "p90", "p95"),
                                            np.percentile(top, [25, 50, 75, 90, 95]), strict=True)),
                                 severe_gt5=sum(top > 5), severe_ge8=sum(top >= 8), severe_gt10=sum(top > 10), empty=len(selected)-len(ok),
                                 wall_s=sum(c["wall_s"] for c in selected),
                                 cpu_s=sum(c["cpu_s"] for c in selected),
                                 evals=sum(c["evals"] for c in selected) if all(c["evals"] is not None for c in selected) else None))
    scan.replay.geom.write_csv(out, "distributions.csv", distribution)
    print(f"{len(cells)} pools, {len(modes)} modes; largest independent RMSD difference {worst:.10g} A")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "verify"))
    for key in ("aurora", "brd9", "out"):
        parser.add_argument(key, type=lambda name: Path(name).resolve())
    for key in ("build", "official", "counter"):
        parser.add_argument("--" + key, type=lambda name: Path(name).resolve())
    args = parser.parse_args()
    if args.action == "run":
        calculate(args)
    analyze(args)
