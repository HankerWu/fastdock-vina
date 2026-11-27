"""Independently recompute assigned RMSDs by Kabsch/SVD and RDKit CalcRMS."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from experiment import IDS, coordinates, reference
from rdkit import Chem
from rdkit.Chem import rdMolAlign


def verify(work):
    with (work / "modes.csv").open() as f:
        expected = {
            (r["key"], int(r["rank"])): float(r["rmsd_A"]) for r in csv.DictReader(f)
        }
    qca, qc, query = reference(work, "1KV2", "B96", "A")
    xyz = query.GetConformer().GetPositions()
    error = []
    for rx, (ccd, _) in IDS.items():
        rca, rc, _ = reference(work, rx, ccd, "A")
        keys = sorted(
            k
            for k, v in qca.items()
            if k in rca
            and np.linalg.norm(v - qc) < 14
            and np.linalg.norm(rca[k] - rc) < 14
        )
        x, y = map(np.stack, ([qca[k] for k in keys], [rca[k] for k in keys]))
        u, _, vt = np.linalg.svd((x - x.mean(0)).T @ (y - y.mean(0)))
        u[:, -1] *= np.sign(np.linalg.det(u @ vt))
        ref = Chem.Mol(query)
        conf = ref.GetConformer()
        for i, v in enumerate((xyz - x.mean(0)) @ u @ vt + y.mean(0)):
            conf.SetAtomPosition(i, v)
        folder = work / ("observed" if rx == "1A9U" else "run/raw")
        for p in folder.glob("*_12000.pdbqt"):
            for idx, block in enumerate(p.read_text().split("ENDMDL")[:-1], 1):
                mol, _, _ = coordinates(block.splitlines())
                rms = rdMolAlign.CalcRMS(
                    mol, ref, maxMatches=10001, symmetrizeConjugatedTerminalGroups=False
                )
                error.append(abs(rms - expected.pop((p.stem, idx))))
    if expected or max(error) > 1e-8:
        raise ValueError((expected, max(error)))
    return {"matched_modes": len(error), "max_rmsd_disagreement_A": max(error)}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    a = p.parse_args()
    (a.work / "verification.json").write_text(
        json.dumps(verify(a.work), indent=2) + "\n"
    )
