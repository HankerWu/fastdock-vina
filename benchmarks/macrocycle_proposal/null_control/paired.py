"""Compare reallocated, locked-null and unrestricted proposals for two registered macrocycles."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path
from statistics import median


def run(a):
    os.sched_setaffinity(0, {a.cpu})
    base = a.repo.resolve() / "benchmarks/macrocycle_proposal"
    sys.path.insert(0, str(base))
    from panel import command, counts, read, write, structural, cases

    sys.path.insert(0, str(base / "external"))
    from external_score import main as evaluate_mer

    pdbs = ("5QCT", "5K0X")
    arms = {"parent": a.parent, "lock": a.lock, "null": a.null}
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    hashes = {arm: digest(p) for arm, p in arms.items()}
    bace = next(c for c in cases(a.inputs, "macro") if c[0] == "5QCT")
    mer = a.holdout / "5K0X"
    receptors = {
        "5QCT": bace,
        "5K0X": (
            "5K0X",
            "K0X",
            mer / "receptor.pdbqt",
            mer / "ligand.pdbqt",
            mer / "box.txt",
        ),
    }
    original = {}
    for d in (a.past_bace, a.past_mer):
        for r in read(d / "runs.csv"):
            if r["pdb"] in pdbs and r["arm"] in {"parent", "lock"}:
                original[r["pdb"], r["arm"], int(r["seed"]), int(r["cap"])] = r
    rows = {}
    times = []
    ranks = []
    cells = list(
        itertools.product(pdbs, (458921, 740317), ((3000, 2000), (12000, 9000)))
    )
    for rep in range(4):
        random.Random(1671 + rep // 2).shuffle(cells)
        for i, (pdb, seed, (cap, low)) in enumerate(cells):
            pair = [
                ("parent", cap),
                ("lock", cap),
                ("null", cap),
                ("lock", low),
                ("null", low),
            ]
            k = i % len(pair)
            pair = pair[k:] + pair[:k]
            if rep % 2:
                pair.reverse()
            d = a.out / pdb
            raw = d / "raw"
            raw.mkdir(exist_ok=True, parents=True)
            for arm, budget in pair:
                key = f"{pdb}_{seed}_{budget}_{arm}"
                pose = raw / (key + ".pdbqt")
                cmd = command(receptors[pdb], seed, budget, arms[arm], pose)
                tic = time.perf_counter_ns()
                log = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
                elapsed = time.perf_counter_ns() - tic
                counters = counts(log)
                sha = digest(pose)
                if arm in {"parent", "lock"}:
                    old = original[pdb, arm, seed, budget]
                    if old.get("sha256", old.get("pool_sha256")) != sha or any(
                        int(old[k]) != v for k, v in counters.items()
                    ):
                        raise ValueError(f"External ledger differs for {key}")
                if rep == 0:
                    (raw / (key + ".log")).write_text(log)
                    rows[key] = dict(
                        key=key,
                        pdb=pdb,
                        seed=seed,
                        cap=budget,
                        arm=arm,
                        modes=pose.read_text().count("ENDMDL"),
                        **counters,
                        joint=sum(
                            v for k, v in counters.items() if k.endswith("_evals")
                        ),
                        sha256=sha,
                    )
                elif rows[key]["sha256"] != sha or any(
                    rows[key][k] != v for k, v in counters.items()
                ):
                    raise ValueError(f"Repeated job changed {key}")
                times.append(dict(rep=rep, key=key, elapsed_ns=elapsed))
    if hashes != {arm: digest(p) for arm, p in arms.items()}:
        raise ValueError("Binary changed during experiment")
    (a.out / "binaries.json").write_text(json.dumps(hashes, indent=2) + "\n")
    write(a.out / "timing.csv", times)
    for pdb in pdbs:
        d = a.out / pdb
        rs = [r for r in rows.values() if r["pdb"] == pdb]
        for r in rs:
            r["median_wall_s"] = (
                median(t["elapsed_ns"] for t in times if t["key"] == r["key"]) * 1e-9
            )
        write(d / "runs.csv", rs)
        if pdb == "5QCT":
            structural(a.inputs, d, rs)
        else:
            evaluate_mer(a.holdout, d)
        structs = read(d / "structures.csv")
        models = read(d / "models.csv")
        for r in structs:
            k = r["key"]
            modes = [m for m in models if m["key"] == k]
            top = float(r.get("top1_A", r.get("top1")))
            best5 = float(r.get("best5_A", r.get("best5")))
            ora = float(r.get("oracle_A", r.get("oracle")))
            bond = float(modes[0].get("cg_distance_A", modes[0].get("closure_bond_A")))
            ranks.append(
                dict(
                    pdb=pdb,
                    key=k,
                    arm=r["arm"],
                    cap=r["cap"],
                    seed=r["seed"],
                    joint=r["joint"],
                    median_wall_s=r["median_wall_s"],
                    top1_A=top,
                    top5_A=best5,
                    oracle_A=ora,
                    closure_A=bond,
                    open_modes=sum(
                        float(m.get("cg_distance_A", m.get("closure_bond_A"))) > 3
                        for m in modes
                    ),
                    modes=len(modes),
                )
            )
    write(a.out / "summary.csv", ranks)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for v in (
        "repo",
        "inputs",
        "holdout",
        "past_bace",
        "past_mer",
        "out",
        "parent",
        "lock",
        "null",
    ):
        p.add_argument(v, type=Path)
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    run(p.parse_args())
