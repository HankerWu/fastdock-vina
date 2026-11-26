"""Crystal-packing sensitivity for monomers with ≤4 Å contacts outside the author assembly."""

import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from pose_assembly import add_mate, mate_transform
from pose_run import HERE, load_meta, measure, write_csv
from run import receptor


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    args = p.parse_args()
    work, out = args.inputs.resolve(), args.output.resolve()
    rows = list(csv.DictReader((HERE / "pose/packing.csv").open()))
    jobs = []
    for r in rows:
        imgs = json.loads(r["close_images"])
        if not imgs or any(i["assemblies"] for i in imgs):
            continue
        pdb = r["pdb"]
        source = work / pdb
        meta = load_meta(work, pdb)
        obs = np.array(meta["mapped"])[0]
        rec0 = source / "receptor.pdbqt"
        dep = source / (pdb + ".pdb")
        tilted = out / "plus_lattice" / pdb
        tilted.mkdir(parents=True, exist_ok=True)
        (tilted / "ligand.pdbqt").write_bytes((source / "ligand.pdbqt").read_bytes())
        assembled = rec0
        for j, im in enumerate(imgs):
            addition = rec0
            if im["chain"] != meta["chain"]:
                tmp = tilted / im["chain"]
                tmp.mkdir(exist_ok=True)
                lines = [
                    l[:21] + "A" + l[22:]
                    for l in dep.read_text().splitlines()
                    if l.startswith("ATOM  ") and l[21] == im["chain"]
                ]
                m, t = mate_transform(dep, im)
                prep = receptor(lines, tmp, pdb, (obs - t) @ m)
                (tmp / "repair.json").write_text(json.dumps(prep) + "\n")
                addition = tmp / "receptor.pdbqt"
            output = tilted / (
                "receptor.pdbqt" if j == len(imgs) - 1 else f"plus{j}.pdbqt"
            )
            add_mate(dep, addition, assembled, im, output)
            assembled = output
        for policy, idir in (("monomer", work), ("plus_lattice", out / "plus_lattice")):
            dest = out / policy / "poses"
            dest.mkdir(parents=True, exist_ok=True)
            for seed in (91571, 91583):
                jobs.append(
                    (
                        policy,
                        (
                            idir,
                            dest,
                            args.binary.resolve(),
                            meta,
                            meta,
                            np.array(meta["mapped"]),
                            0.0,
                            0,
                            0,
                            seed,
                            "armijo",
                            120000,
                        ),
                    )
                )

    def one(job):
        policy, args = job
        cell, modes = measure(args)
        return dict(policy=policy, **cell), [dict(policy=policy, **m) for m in modes]

    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(one, jobs))
    write_csv(out / "cells.csv", [v for v, _ in results])
    write_csv(out / "modes.csv", [m for _, modes in results for m in modes])


if __name__ == "__main__":
    main()
