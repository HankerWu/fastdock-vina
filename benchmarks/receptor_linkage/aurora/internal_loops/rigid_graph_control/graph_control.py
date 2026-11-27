"""Isolate the graph interface while retaining the Pro release and counters."""

import argparse
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "rigid_coupled"))
import study
from study import scan, audit


def run(a):
    raw = a.out / "raw"
    raw.mkdir(parents=True, exist_ok=False)
    cases = scan.config(a.aurora, a.brd9, a.out)
    manifest = json.loads(
        (HERE.parent / "mutation_scan/results/inputs.json").read_text()
    )
    rows, counts, commands = [], [], []
    for drug, c in cases.items():
        for key in ("ligand", "bonds", "receptor"):
            assert scan.replay.digest(c[key]) == manifest[drug][key]
        for i, seed in enumerate(study.SEEDS):
            for arm in (("pro", "pivot") if i == 0 else ("pivot", "pro")):
                exe = a.build / ("base" if arm == "pro" else "pivot") / "vina"
                stem = f"{drug}_{seed}_{arm}"
                pool = raw / (stem + ".pdbqt")
                log = pool.with_suffix(".log")
                cmd = scan.command(exe, c, seed, pool)
                k = cmd.index("--ligand_bonds")
                del cmd[k : k + 2]
                cost = scan.replay.invoke(cmd, log)
                commands.append(cmd)
                row = dict(drug=drug, seed=seed, arm=arm)
                rows.append(
                    dict(
                        **row,
                        **cost,
                        binary_sha256=scan.replay.digest(exe),
                        pool_sha256=scan.replay.digest(pool),
                    )
                )
                counts.append(dict(**row, **study.totals(log)))
    scan.replay.geom.write_csv(a.out, "jobs.csv", rows)
    scan.replay.geom.write_csv(a.out, "counts.csv", counts)
    study.dump(a.out / "commands.json", commands)
    study.dump(a.out / "inputs.json", manifest)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("run", "verify"))
    for k in ("aurora", "brd9", "out"):
        p.add_argument(k, type=lambda s: Path(s).resolve())
    p.add_argument("--build", type=lambda s: Path(s).resolve())
    a = p.parse_args()
    if a.action == "run":
        run(a)
    # Reuse the main panel's independent coordinate readers and distributions, for these arms.
    study.ARMS = {k: study.ARMS[k] for k in ("pro", "pivot")}
    study.analyze(a)
