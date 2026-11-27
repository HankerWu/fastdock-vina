"""Measure shared input-preparation cost from local source files and check exact inputs."""

import argparse
import hashlib
import os
import shutil
import tempfile
import time
from pathlib import Path
from inputs import CASES, prepare, reference


def measure(work):
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        dst = Path(tmp)
        (dst / "raw").mkdir()
        for pdb, (_, ccd, _) in CASES.items():
            for name in (pdb + ".pdb", ccd + ".cif", "CL.cif", "GOL.cif", "SO4.cif"):
                shutil.copyfile(work / "raw" / name, dst / "raw" / name)
        for pdb in CASES:
            t = time.perf_counter()
            prepare(dst, pdb)
            elapsed = time.perf_counter() - t
            hashes = {}
            for name in ("ligand", "receptor"):
                out = (dst / pdb / (name + ".pdbqt")).read_bytes()
                hashes[name + "_sha256"] = hashlib.sha256(out).hexdigest()
                if out != (work / pdb / (name + ".pdbqt")).read_bytes():
                    raise ValueError((pdb, name))
            rows.append(dict(pdb=pdb, seconds=elapsed, **hashes))
    reference.write(work / "preparation_cost.csv", rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument("--cpu", type=int, default=5)
    args = ap.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    measure(args.work)
