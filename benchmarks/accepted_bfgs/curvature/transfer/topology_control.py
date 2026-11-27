"""Replay just the changed input graph in two solvers; VH032 is the exact null."""
import argparse
import csv
import json
import resource
import sys
import subprocess
import time
from pathlib import Path

import numpy as np

from evaluate import PDBS, frame, sha, write_csv

BENCH = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(BENCH / "receptor_linkage/aurora"))
from geometry import mode_coords


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def run(work, out, binaries):
    out.mkdir(parents=True, exist_ok=True)
    table = Path(__file__).resolve().parent / "results"
    previous = {(r["query"], r["receptor"], r["seed"], r["arm"]): r
                for r in read(table / "cells.csv")}
    meta = {p:json.loads((work / p / "metadata.json").read_text()) for p in PDBS}
    rows, modes = [], []
    for cell in read(table / "cost.csv"):
        q, r, seed, arm = cell["key"].split("_", 3)
        if arm not in binaries:
            continue
        args = json.loads(cell["arguments"])
        original = previous[q,r,seed,arm]
        inp = work / args[args.index("--out") + 1]
        if sha(inp) != original["output_sha256"]:
            raise ValueError(f"Unexpected original pool: {inp}")
        dest = out / (cell["key"] + ".pdbqt")
        args[args.index("--out") + 1] = dest
        # Bonds use only the query's chemical graph and the prepared index map.
        bonds = table.parent / "bond_graphs" / f"{q}.bonds"
        cpu0 = resource.getrusage(resource.RUSAGE_CHILDREN)
        t0 = time.perf_counter()
        with dest.with_suffix(".log").open("w") as lg:
            subprocess.run([binaries[arm], *args, "--ligand_bonds", bonds],
                           cwd=work, stdout=lg, stderr=subprocess.STDOUT, check=True)
        wall = time.perf_counter() - t0
        cpu = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_s = cpu.ru_utime + cpu.ru_stime - cpu0.ru_utime - cpu0.ru_stime
        ref, _, _ = frame(meta[q],meta[r])
        coors = mode_coords(dest, meta[q]["serial"])
        rms = np.sqrt(np.square(coors[:,None] - ref).sum(-1).mean(-1).min(-1))
        same = inp.read_bytes() == dest.read_bytes()
        if q == "4W9H" and not same:
            raise ValueError("Chemically identical input graph changed output coordinates")
        for k, val in enumerate(rms, 1):
            modes.append(dict(query=q,receiver=r,seed=seed,arm=arm,rank=k,rmsd_A=val))
        rows.append(dict(query=q,receiver=r,seed=seed,arm=arm,blind_sha256=sha(inp),
                         graph_locked_sha256=sha(dest),identical=same,modes=len(rms),
                         old_top1_A=original["top1"],new_top1_A=rms[0],
                         old_best5_A=original["best5"],new_best5_A=min(rms[:5]),
                         old_oracle_A=original["oracle"],new_oracle_A=min(rms),
                         old_regret_A=original["regret"],new_regret_A=rms[0]-min(rms),
                         wall_s=wall,cpu_s=cpu_s))
    if len(rows)!=8:
        raise ValueError("Incomplete normal-budget panel")
    write_csv(out / "cells.csv",rows)
    write_csv(out / "modes.csv",modes)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work","out","parent","fixed"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    run(a.work.resolve(),a.out.resolve(),{"parent":a.parent.resolve(),"fixed":a.fixed.resolve()})
