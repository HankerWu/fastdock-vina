"""Independent S6K1 pocket audit and cross-check against the prior production rigid-window evaluator."""

import argparse
import csv
import io
import json
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


def audit(pdb):
    lines = pdb.read_text().splitlines()

    def xyz(line):
        return [float(line[j : j + 8]) for j in [30, 38, 46]]

    het = [
        l
        for l in lines
        if l.startswith("HETATM")
        and l[16] in " A"
        and l[76:78].strip() not in ("H", "D")
    ]
    stu = np.array([xyz(l) for l in het if l[17:20] == "STU" and l[21] == "A"])
    atoms = [
        l
        for l in lines
        if l.startswith("ATOM") and l[16] in " A" and l[76:78].strip() not in ("H", "D")
    ]

    def nearest(records):
        return (
            round(float(cdist(stu, [xyz(l) for l in records]).min()), 3)
            if records
            else None
        )

    missing = sorted(
        int(m[1])
        for line in lines
        if (m := re.match(r"REMARK 465\s+[A-Z]{3}\s+A\s+(\d+)\s*$", line))
    )
    stretches = []
    for residue in missing:
        if stretches and stretches[-1][-1] == residue - 1:
            stretches[-1].append(residue)
        else:
            stretches.append([residue])
    gaps = []
    for s in stretches:
        flank = [
            l for l in atoms if l[21] == "A" and int(l[22:26]) in [s[0] - 1, s[-1] + 1]
        ]
        gaps.append(dict(start=s[0], end=s[-1], flank_min_A=nearest(flank)))
    return dict(
        chain="A",
        donor="STU",
        missing_residues=gaps,
        nearest_crystal_water_A=nearest([l for l in het if l[17:20] == "HOH"]),
        nearest_other_hetero_A=nearest(
            [
                l
                for l in het
                if l[17:20] != "HOH" and not (l[17:20] == "STU" and l[21] == "A")
            ]
        ),
        nearest_other_chain_protein_A=nearest([l for l in atoms if l[21] != "A"]),
    )


def compare(a):
    p = pd.read_csv(a.s6k / "sensitivity.csv").query("nodes==8192")
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        dest = Path(tmp) / "top.pdbqt"
        for r in p.itertuples(index=False):
            raw = (a.s6k / "pools" / f"3A60_{r.ligand}_{r.seed}_pro.pdbqt").read_text()
            dest.write_text(raw.split("ENDMDL", 1)[0] + "ENDMDL\n")
            args = [
                a.binary,
                "--score_only",
                "--receptor",
                a.s6k / "prepared/3A60.pdbqt",
                "--ligand",
                a.s6k / "prepared" / (str(r.ligand) + ".pdbqt"),
                "--rigid_window_pool",
                dest,
                "--rigid_window_samples",
                "4096",
                "--cpu",
                "1",
                "--verbosity",
                "0",
            ]
            proc = subprocess.run(args, text=True, capture_output=True)
            prefix = "model,inter,intra,slope,"
            values = dict(ligand=r.ligand, seed=r.seed, sobol_rigid=r.rigid)
            if proc.returncode:
                if "closure or dummy atoms" not in (proc.stdout + proc.stderr).lower():
                    raise RuntimeError(proc.stderr or proc.stdout)
                values["status"] = "production excludes macrocycle"
            else:
                table = csv.DictReader(
                    io.StringIO(prefix + proc.stdout.split(prefix, 1)[1])
                )
                row = next(table)
                values.update(
                    status="ok",
                    production_rigid=float(row["window_delta"]),
                    center_error=float(row["inter"]) - r.inter,
                )
            results.append(values)
    pd.DataFrame(results).to_csv(a.s6k / "native_rigid_audit.csv", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("s6k", type=Path)
    p.add_argument("binary", type=Path)
    a = p.parse_args()
    (a.s6k / "structure_audit.json").write_text(
        json.dumps(audit(a.s6k / "raw/3A60.pdb"), indent=2) + "\n"
    )
    compare(a)
