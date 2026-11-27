"""Check the archived unrounded stage witnesses against independent crystal maps."""
import argparse
import csv
import json
import tempfile
from pathlib import Path

import numpy as np
from rdkit.Chem import rdMolAlign

import study


def verify(aurora, brd9, observed):
    with tempfile.TemporaryDirectory() as d:
        systems = study.audit.inputs(aurora, brd9, Path(d))
    rows = list(csv.DictReader((observed / "stages.csv").open()))
    first = list(csv.DictReader((observed / "splits.csv").open()))
    maximum = {"stage": 0., "fork": 0.}
    with np.load(observed / "witnesses.npz", allow_pickle=False) as w:
        for group, records in (("stage", rows), ("fork", first)):
            for row in records:
                drug, seed, arm = (row[k] for k in ("drug", "seed", "arm"))
                mol, serial, ref = systems[drug][-3:]
                reorder = [w[drug + "_serials"].tolist().index(i) for i in serial]
                if group == "stage":
                    keys = {f"{drug}_{seed}_{arm}_{row['stage']}": "nearest_A"}
                else:
                    prefix = f"{drug}_{seed}_{row['trajectory']}_{arm}_"
                    keys = {prefix + k: col for k, col in (("parent", "parent_A"), ("proposal", "start_A"), ("return", "return_A"))}
                for key, column in keys.items():
                    for i, position in enumerate(w[key][reorder]):
                        mol.GetConformer().SetAtomPosition(i, position)
                    rms = rdMolAlign.CalcRMS(mol, ref, symmetrizeConjugatedTerminalGroups=False)
                    maximum[group] = max(maximum[group], abs(rms - float(row[column])))
    assert max(maximum.values()) < 1e-5
    print(json.dumps(dict(stages=len(rows), fork_states=len(first)*3, max_abs_A=maximum), indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("aurora", "brd9", "observed"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    verify(a.aurora, a.brd9, a.observed)
