"""Replay both 3P6E cognate merge streams and verify all ordered inputs agree across arms."""

import csv
import gzip
import json
from pathlib import Path
import sys

import numpy as np

from transfer_stages import replay


def stream(path):
    with gzip.open(path, "rt") as f:
        lines = list(f)
    serial = [
        int(i)
        for i in next(line for line in lines if line.startswith("atoms")).split()[1:]
    ]
    return serial, [line.split("\t")[0] for line in lines if line.startswith("merge ")]


def read(raw, inputs):
    meta = json.loads((inputs / "3P6E/preparation.json").read_text())
    reference = np.asarray(meta["mapped"])
    for seed in (684173, 961751):
        base = f"3P6E_3P6E_{seed}_"
        old = stream(raw / (base + "parent.trace.gz"))
        if old != stream(raw / (base + "candidate.trace.gz")):
            raise ValueError("Unmatched ordered merge inputs")
        atoms, proposals = old
        order = [atoms.index(i) for i in meta["serial"]]
        xyz = {}
        previous, fixed = [], []
        for step, line in enumerate(proposals, 1):
            _, radius, cap, energy, points = line.split(" ", 4)
            pos = np.fromstring(points, sep=" ").reshape(-1, 3)
            key = (float(energy), *pos[0])
            xyz[key] = pos
            previous = replay.update(
                previous, key, float(radius), int(cap), xyz, "nearest"
            )
            fixed = replay.update(fixed, key, float(radius), int(cap), xyz, "all")
            rmsd = np.sqrt(
                np.square(reference - pos[order]).sum(axis=2).mean(axis=1)
            ).min()
            yield (
                seed,
                step,
                float(energy),
                rmsd,
                int(key in previous),
                int(key in fixed),
                len(previous),
                len(fixed),
            )


if __name__ == "__main__":
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(
        (
            "seed",
            "arrival",
            "energy",
            "rmsd_A",
            "parent_kept",
            "fixed_kept",
            "parent_size",
            "fixed_size",
        )
    )
    writer.writerows(read(Path(sys.argv[1]), Path(sys.argv[2])))
