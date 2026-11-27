"""Preserve the tags, scalar energy and receptor-frame error for every trace row."""

import argparse
import csv
from pathlib import Path

import numpy as np

from branches import observer, study


def archive(aurora, brd9, trace_dir, destination):
    destination.mkdir(parents=True, exist_ok=False)
    systems = study.scan.config(aurora, brd9, trace_dir)
    with (trace_dir / "observation_cost.csv").open() as f:
        for row in csv.DictReader(f):
            stem = "_".join(row[k] for k in ("drug", "seed", "arm"))
            path = trace_dir / (stem + ".bin")
            assert study.scan.replay.digest(path) == row["raw_sha256"]
            serials, tr = observer.read_trace(path)
            system = systems[row["drug"]]
            idx = [list(system["serials"]).index(s) for s in serials]
            rms = observer.distances(tr, system["refs"][:, idx])
            fields = {key: tr[key] for key in ("traj", "step", "phase", "event", "local", "energy")}
            np.savez_compressed(destination / (stem + ".npz"), rmsd_A=rms, **fields)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("aurora", "brd9", "trace_dir", "destination"):
        p.add_argument(key, type=Path)
    a = p.parse_args()
    archive(a.aurora, a.brd9, a.trace_dir, a.destination)
