"""Reproduce the registered three-arm comparison; parent and fixed get identical inputs."""
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("inputs", "output", "upstream", "upstream_count", "parent", "fixed"):
        parser.add_argument(arg, type=Path)
    parser.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    parser.add_argument("--diagnostic-boxes", type=Path, help="Replay only the two observed marker-face stress boxes")
    a = parser.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    engines = {arm: getattr(a, arm).resolve() for arm in ("upstream", "parent", "fixed")}
    raw = a.output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    rows = []
    ids = ("3VHA", "6QBH") if a.diagnostic_boxes else ("3VHA", "6QBH", "8H6P")
    jobs = list(itertools.product(ids, (191009, 458921), (3000, 12000)))
    random.Random(51829).shuffle(jobs)
    for i, (pdb, seed, cap) in enumerate(jobs):
        folder = a.inputs / pdb
        center = json.loads((folder / "preparation.json").read_text())["center"]
        box = a.output / f"{pdb}.box"
        if a.diagnostic_boxes:
            box.write_bytes((a.diagnostic_boxes / (pdb + ".box")).read_bytes())
        else:
            box.write_text("".join(f"center_{d} = {v}\nsize_{d} = 26\n" for d, v in zip("xyz", center)))
        case = (pdb, "", folder / "receptor.pdbqt", folder / "ligand.pdbqt", box)
        labels = list(engines)
        for arm in labels[i % 3:] + labels[:i % 3]:
            key = f"{pdb}_{seed}_{cap}_{arm}"
            pose, log = raw / (key + ".pdbqt"), raw / (key + ".log")
            wall, cpu, rss = timed(command(case, seed, cap, engines[arm], pose, arm != "upstream"), log)
            if arm == "upstream":
                replay = raw / (key + ".repeat")
                log = raw / (key + ".counts")
                timed(command(case, seed, cap, a.upstream_count, replay, False), log)
                assert pose.read_bytes() == replay.read_bytes()
                replay.unlink()
            rows.append(dict(key=key, pdb=pdb, seed=seed, cap=cap, arm=arm, **counts(log.read_text()),
                             wall_s=wall / 1e9, cpu_s=cpu, max_rss_KiB=rss, modes=pose.read_text().count("ENDMDL"),
                             sha256=sha(pose)))
    write(a.output / "runs.csv", rows)
    evaluate(a.inputs, a.output)
    engines["upstream_count"] = a.upstream_count
    (a.output / "binaries.json").write_text(json.dumps({k: sha(v) for k, v in engines.items()}, indent=2) + "\n")


if __name__ == "__main__":
    main()
