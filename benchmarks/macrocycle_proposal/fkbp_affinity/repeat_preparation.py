"""Repeat cached source-file preparation, verify input hashes and CCD stereo."""

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
for arg in ("repo", "work", "out"):
    p.add_argument(arg, type=Path)
a = p.parse_args()
sys.path.insert(0, str(a.repo.resolve() / "benchmarks/bond_distance"))
from chem_graph import check


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


rows = []
with tempfile.TemporaryDirectory() as tmp:
    target = Path(tmp)
    shutil.copytree(a.work / "raw", target / "raw")
    for pdb, component in (("1FKJ", "FK5"), ("6VCU", "R27"), ("1FKF", "FK5")):
        start = time.perf_counter_ns()
        subprocess.run(
            [
                sys.executable,
                Path(__file__).with_name("prepare_assay.py"),
                a.repo,
                target,
                pdb,
            ],
            check=True,
        )
        elapsed = time.perf_counter_ns() - start
        g, mapping, ends, centers, alkenes, mismatch = check(
            target / pdb / "ligand.pdbqt", target / "raw" / (component + ".cif")
        )
        files = ["ligand.pdbqt"]
        if pdb != "6VCU":
            files.append("receptor.pdbqt")
        for fn in files:
            if digest(a.work / pdb / fn) != digest(target / pdb / fn):
                raise RuntimeError(f"Input changed: {pdb}/{fn}")
        rows.append(
            dict(
                pdb=pdb,
                ccd=component,
                pdb_sha256=digest(target / "raw" / (pdb + ".pdb")),
                ccd_sha256=digest(target / "raw" / (component + ".cif")),
                ligand_sha256=digest(target / pdb / "ligand.pdbqt"),
                receptor_sha256=digest(target / pdb / "receptor.pdbqt")
                if pdb != "6VCU"
                else "",
                preparation_s=elapsed / 1e9,
                heavy_atoms=len(mapping),
                formal_charge=sum(x.GetFormalCharge() for x in g.GetAtoms()),
                stereocenters=centers,
                specified_double_bonds=alkenes,
                mismatches=len(mismatch),
                closure_atoms=" ".join(sorted(ends)),
            )
        )
with (a.out / "input_verification.csv").open("w") as f:
    w = csv.DictWriter(f, rows[0], lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
