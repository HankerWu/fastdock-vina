"""Check geometry independently using Meeko SMILES serials and PDB-derived graphs."""

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from run_transfer import BENCH, GROUPS

sys.path.insert(0, str(BENCH / "paired_crossdocking"))
import independent_cross as independent


def verify(work):
    independent.COMP = {s: c for entries in GROUPS.values() for s, c in entries.items()}
    inputs = {}
    for s, comp in independent.COMP.items():
        ca, _, refs, serial, other = independent.complex_data(work, s)
        ls = (work / "raw" / f"{s}.pdb").read_text().splitlines()
        points = [
            independent.xyz(l)
            for l in ls
            if l.startswith("HETATM")
            and l[21] == "A"
            and l[17:20].strip() == comp
            and l[76:78].strip() not in ("H", "D")
        ]
        inputs[s] = ca, np.mean(points, axis=0), refs, serial, other
    geometry = pd.read_csv(work / "geometry.csv").set_index(
        ["query", "receiver", "cap", "seed", "engine", "model"]
    )
    count = 0
    worst = 0.0
    for pool in sorted((work / "runs").glob("*.pdbqt")):
        q, r, cap, seed, engine = pool.stem.split("_")
        ref, _, n, _ = independent.fitted(inputs[q], inputs[r])
        for b in pool.read_text().split("ENDMDL"):
            if not b.strip():
                continue
            lines = b.splitlines()
            rank = int(next(l for l in lines if l.startswith("MODEL")).split()[1])
            coords = {
                int(l[6:11]): independent.xyz(l)
                for l in lines
                if l.startswith(("ATOM", "HETATM"))
            }
            pos = np.array([coords[s] for s in inputs[q][3]])
            rms = min(np.linalg.norm(x - pos) / np.sqrt(len(pos)) for x in ref)
            row = geometry.loc[(q, r, int(cap), int(seed), engine, rank)]
            if n != row.n_ca:
                raise ValueError(f"Inconsistent receptor fit: {q} {r}")
            worst = max(worst, abs(rms - row.rmsd))
            count += 1
    if count != len(geometry) or worst > 1.0e-7:
        raise ValueError((count, len(geometry), worst))
    print(
        f"{count} heavy-atom RMSDs verified, maximum absolute discrepancy {worst:.3g} Å."
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    verify(p.parse_args().work)
