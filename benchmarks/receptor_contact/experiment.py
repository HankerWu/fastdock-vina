"""Frozen receptor-directed pilot on the BRD4 first bromodomain crystallographic panel."""

import argparse
import csv
import itertools
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmarks/mapping_ensemble"))
from dock_matrix import run_one

LIGANDS = ("4O76", "4O77", "4OGI", "4OGJ")
RECEPTORS = ("3MXF", "4OGI", "4OGJ")
SEEDS = (314159, 271828, 161803, 141421)
ARMS = ("native", "multi", "contact")
DISTANCE = 2.9


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def prepare(data, out):
    out.mkdir(exist_ok=True, parents=True)
    inventory = []
    for rec in RECEPTORS:
        lines = (data / f"{rec}_receptor.pdbqt").read_text().splitlines()
        res = {
            l[12:16].strip(): np.array([float(l[k : k + 8]) for k in (30, 38, 46)])
            for l in lines
            if l.startswith("ATOM  ") and l[17:20] == "ASN" and int(l[22:26]) == 140
        }
        nd = res["ND2"]
        for lig in LIGANDS:
            atoms = (data / f"{lig}_ligand.pdbqt").read_text().splitlines()
            acceptors = [
                int(l[6:11])
                for l in atoms
                if l.startswith(("ATOM  ", "HETATM"))
                and l.split()[-1] in ("NA", "OA", "SA")
            ]
            if not acceptors:
                raise ValueError(f"No acceptor in {lig}")
            blocks = []
            for h in ("HD21", "HD22"):
                v = res[h] - nd
                if not 0.85 < np.linalg.norm(v) < 1.15:
                    raise ValueError("Invalid amide geometry")
                target = nd + DISTANCE * v / np.linalg.norm(v)
                for serial in acceptors:
                    blocks.append(f"{serial} " + " ".join(f"{x:.8f}" for x in target))
                    inventory.append(
                        {
                            "receptor": rec,
                            "ligand": lig,
                            "hydrogen": h,
                            "serial": serial,
                            "x": target[0],
                            "y": target[1],
                            "z": target[2],
                        }
                    )
            (out / f"{rec}_{lig}.contact").write_text("\nOR\n".join(blocks) + "\n")
    write(out / "sites.csv", inventory)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("prepared", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", type=int, default=4)
    a = p.parse_args()
    prepare(a.prepared, a.work)
    with ThreadPoolExecutor(a.jobs) as pool:
        list(
            pool.map(
                lambda c: run_one(a.binary.resolve(), a.prepared, a.work, c),
                itertools.product(RECEPTORS, LIGANDS, SEEDS, ("contact",)),
            )
        )
