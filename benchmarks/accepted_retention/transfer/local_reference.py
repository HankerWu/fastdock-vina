"""Diagnostic relaxation of deposited and torsion-projected geometries after the blind dock."""

import argparse
import csv
import json
import sys
from pathlib import Path
from inputs import ROOT, CASES, reference

sys.path.insert(0, str(ROOT / "paired_crossdocking"))
from crystal_start import run

sys.path.insert(0, str(ROOT / "pair_gradient"))
from local_basins import common_score

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument("parent", type=Path)
    a = ap.parse_args()
    m = {
        p: json.loads((a.work / p / "preparation.json").read_text())
        for p in sorted(CASES)
    }
    (a.work / "metadata.json").write_text(json.dumps(m))
    run(a.work.resolve(), a.parent.resolve())

    with (a.work / "crystal_starts.csv").open() as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        q, r, kind = row["source"], row["receiver"], row["arm"]
        common = [
            "--receptor",
            a.work / r / "receptor.pdbqt",
            "--seed",
            "43103",
            "--cpu",
            "1",
        ]
        for axis, x in zip("xyz", m[r]["center"]):
            common.extend(["--center_" + axis, str(x), "--size_" + axis, "24"])
        ligand = a.work / q / "ligand.pdbqt"
        raw = a.work / "crystal_diagnostic"
        row["common_inter_plus_intra"] = common_score(
            a.parent, common, ligand, raw / f"{q}-{r}-{kind}.pdbqt"
        )
        row["common_initial_inter_plus_intra"] = common_score(
            a.parent, common, ligand, raw / f"{q}-{r}-{kind}-input.pdbqt"
        )
    reference.write(a.work / "crystal_starts.csv", rows)
