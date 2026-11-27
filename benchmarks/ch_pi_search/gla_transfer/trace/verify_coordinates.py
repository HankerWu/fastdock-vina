"""Recompute RMSD from every losslessly archived trajectory record by atom serial."""

import argparse
import csv
import hashlib
import io
import json
import lzma
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from verify import geometry


def verify(work, tables, archive):
    refs_by_run = defaultdict(dict)
    with lzma.open(tables / "events.csv.xz", "rt") as stream:
        for row in csv.DictReader(stream):
            key = (row["stage"], int(row["trajectory"]), int(row["index"]))
            run = refs_by_run[row["run"]]
            assert key not in run
            run[key] = float(row["rmsd_A"])
    with (tables / "replay.csv").open() as stream:
        hashes = {r["tag"]: r for r in csv.DictReader(stream)}
    maximum = 0.0
    matched = 0
    with tarfile.open(archive, "r:xz") as bundle:
        payloads = {
            m.name.rsplit("/", 1)[-1]: bundle.extractfile(m).read()
            for m in bundle
            if m.isfile()
        }
        for tag, expected in refs_by_run.items():
            raw = payloads[tag + ".trace"]
            assert hashlib.sha256(raw).hexdigest() == hashes[tag]["trace_sha256"]
            assert (
                hashlib.sha256(payloads[tag + ".pdbqt"]).hexdigest()
                == hashes[tag]["pool_sha256"]
            )
            panel, query, receptor, seed, cap, search = tag.split("_")
            refs, order = geometry(work, panel, query, receptor)
            lines = raw.splitlines()
            (header,) = [line for line in lines if line.startswith(b"atoms ")]
            positions = [list(map(int, header.split()[1:])).index(i) for i in order]
            rows = [
                l.decode().split(maxsplit=4)
                for l in lines
                if not l.startswith(b"atoms ")
            ]
            xyz = np.loadtxt(io.StringIO("\n".join(r[4] for r in rows))).reshape(
                len(rows), len(order), 3
            )
            rmsd = np.sqrt(
                np.square(xyz[:, None, positions, :] - refs[None, :, :, :])
                .sum(-1)
                .mean(-1)
            ).min(axis=1)
            keys = [(r[0], int(r[1]), int(r[2])) for r in rows]
            assert len(keys) == len(set(keys)) == len(expected)
            old = np.array([expected[key] for key in keys])
            maximum = max(maximum, float(np.abs(old - rmsd).max()))
            matched += int(((old <= 2) == (rmsd <= 2)).sum())
    total = sum(map(len, refs_by_run.values()))
    assert set(hashes) == set(refs_by_run) and matched == total and maximum < 1e-12
    return dict(
        examined=total,
        classification_matches=matched,
        maximum_rmsd_difference_A=maximum,
        verified_trajectory_and_released_pool_sha256=len(hashes),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("tables", type=Path)
    p.add_argument("archive", type=Path)
    a = p.parse_args()
    print(json.dumps(verify(a.work, a.tables, a.archive), indent=2))
