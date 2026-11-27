"""Run all registered paired native searches and score both deposited-coordinate maps."""
import argparse
import itertools
import json
import os
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from benchmarks.concerted_ring.run import command, counts, evaluate, sha, timed, write

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "output", "upstream", "parent", "fixed"):
        p.add_argument(name, type=Path)
    p.add_argument("--higher-budget", action="store_true")
    a = p.parse_args()
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    engines = {k: getattr(a, k).resolve() for k in ("upstream", "parent", "fixed")}
    raw = a.output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    budgets = (48000,) if a.higher_budget else (3000, 12000)
    jobs = list(itertools.product(("3VHA", "4EGK", "2J9M"), (191009, 458921), budgets))
    random.Random(51829).shuffle(jobs)
    rows = []
    for i, (pdb, seed, cap) in enumerate(jobs):
        directory = a.inputs / pdb
        if json.loads((directory / "screening.json").read_text())["status"] != "ready":
            raise ValueError("Ineligible structure " + pdb)
        center = json.loads((directory / "preparation.json").read_text())["center"]
        box = a.output / (pdb + ".box")
        box.write_text("".join(f"center_{d} = {v}\nsize_{d} = 26\n" for d,v in zip("xyz",center)))
        case = (pdb, "", directory / "receptor.pdbqt", directory / "ligand.pdbqt", box)
        arms = list(engines)
        for arm in arms[i%3:] + arms[:i%3]:
            key = f"{pdb}_{seed}_{cap}_{arm}"
            pose, log = raw / (key + ".pdbqt"), raw / (key + ".log")
            elapsed, cpu, rss = timed(command(case, seed, cap, engines[arm], pose, arm != "upstream"), log)
            ctr = counts(log.read_text()) if arm != "upstream" else dict.fromkeys(
                ("hunt_calls","hunt_evals","second_calls","second_evals","explicit_calls","explicit_evals"), "")
            rows.append(dict(key=key, pdb=pdb, seed=seed, cap=cap, arm=arm, **ctr,
                             wall_s=elapsed/1e9, cpu_s=cpu, max_rss_KiB=rss,
                             modes=pose.read_text().count("ENDMDL"), sha256=sha(pose)))
    write(a.output / "runs.csv", rows)
    evaluate(a.inputs, a.output)
    (a.output / "binaries.json").write_text(json.dumps({k:sha(v) for k,v in engines.items()}, indent=2)+"\n")
