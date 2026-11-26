"""Independently verify archived basin energies on the original and reloaded starting graphs."""

import argparse
import csv
import json
import math
from pathlib import Path
import re
import subprocess


def score(vina, case, center, ligand, pose):
    cmd = [
        str(vina),
        "--receptor",
        str(case / "receptor.pdbqt"),
        "--ligand",
        str(ligand),
        "--pose_coordinates",
        str(pose),
        "--score_only",
        "--unbound_energy",
        "0",
        "--seed",
        "63419",
        "--cpu",
        "1",
    ]
    for axis, value in zip("xyz", center):
        cmd += ["--center_" + axis, str(value), "--size_" + axis, "24"]
    stdout = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    return [
        float(re.search(key + r"\s*:\s*([-+.\d]+)", stdout)[1])
        for key in ("Final Intermolecular Energy", "Final Total Internal Energy")
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixtures", type=Path)
    parser.add_argument(
        "poses", type=Path, help="extracted starts_and_relaxations.tar.xz"
    )
    parser.add_argument("vina", type=Path)
    args = parser.parse_args()
    max_rms = 0
    rows = list(csv.DictReader(Path(__file__).with_name("energies.csv").open()))
    for row in rows:
        case = args.fixtures / row["pdb"]
        meta = json.loads((case / "metadata.json").read_text())
        start = args.poses / row["pdb"] / (row["start"] + ".pdbqt")
        pose = (
            start.with_name(start.stem + "-relaxed.pdbqt")
            if row["stage"] == "relaxed"
            else start
        )
        xyz = {
            int(line[6:11]): [float(line[k : k + 8]) for k in (30, 38, 46)]
            for line in pose.read_text().splitlines()
            if line.startswith(("ATOM  ", "HETATM"))
        }
        rmsd = min(
            math.sqrt(
                math.fsum(
                    math.dist(xyz[s], p) ** 2 for s, p in zip(meta["serial"], mapping)
                )
                / len(mapping)
            )
            for mapping in meta["mapped"]
        )
        max_rms = max(max_rms, abs(rmsd - float(row["rmsd"])))
        for ligand in (case / "ligand.pdbqt", start):
            if score(args.vina.resolve(), case, meta["center"], ligand, pose) != [
                float(row["inter"]),
                float(row["intra"]),
            ]:
                raise ValueError(f"{pose}: energy differs at serialized precision")
    if max_rms > 1e-8:
        raise ValueError("Coordinate audit mismatch")
    print(
        f"{len(rows)} coordinate records and {2 * len(rows)} fixed scores verified; max RMSD difference {max_rms:g} Å"
    )


if __name__ == "__main__":
    main()
