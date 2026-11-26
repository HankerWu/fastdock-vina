"""Cross-check transfer docking RMSDs and native runs against independent public structures."""

import argparse
import csv
import hashlib
import sys
from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmarks/cutoff_switch"))
from verify_cutoff_poses import verify


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def crosscheck(panel, data, runs):
    modes = defaultdict(list)
    for r in table(runs / "modes.csv"):
        modes[(r["ligand"], r["receptor"], r["seed"], r["arm"])].append(r)
    cells = table(runs / "cells.csv")
    archived = ROOT / "benchmarks/cutoff_switch"
    if panel == "brd4":
        checksums = dict(
            s.split()[::-1]
            for s in (archived / "brd4/poses.sha256").read_text().splitlines()
        )
    else:
        checksums = {
            f"{r['ligand']}_{r['receptor']}_{r['seed']}_{r['arm']}.pdbqt": r["sha256"]
            for r in table(archived / "mixed_t4l.csv")
        }
    n_native = 0
    with TemporaryDirectory() as temporary:
        temp = Path(temporary)
        (temp / "raw").symlink_to(data.resolve() / "raw", target_is_directory=True)
        (temp / "runs").mkdir()
        for cell in cells:
            lig, rec, seed, arm = (
                cell[k] for k in ("ligand", "receptor", "seed", "arm")
            )
            name = f"{lig}_{rec}_{seed}_{arm}.pdbqt"
            source = (runs / f"{rec}_{lig}_{seed}_{arm}.pdbqt").resolve()
            (temp / "runs" / name).symlink_to(source)
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            if digest != cell["sha256"]:
                raise ValueError(f"Output digest mismatch for {source}")
            if arm == "native":
                if checksums[name] != digest:
                    raise ValueError(f"Historical native search changed for {name}")
                n_native += 1
            cell["pool_rmsds"] = ",".join(
                r["rmsd"] for r in modes[(lig, rec, seed, arm)]
            )
        with (temp / "cells.csv").open("w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=cells[0], lineterminator="\n")
            writer.writeheader()
            writer.writerows(cells)
        verify(temp)
    print(f"{n_native} native outputs match prior run SHA-256")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("panel", choices=("brd4", "t4l"))
    p.add_argument(
        "data", type=Path, help="Panel directory with original RCSB structures in raw/"
    )
    p.add_argument("runs", type=Path)
    a = p.parse_args()
    crosscheck(a.panel, a.data, a.runs)
