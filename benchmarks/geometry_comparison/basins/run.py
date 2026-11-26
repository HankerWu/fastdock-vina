"""Re-score projected experimental and emitted starts on one prepared ligand graph."""

import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "paired_crossdocking"))
from crystal_start import read_pose, rms
from projection import project


def write_coordinates(template, other, path):
    """Copy only the coordinate fields and verify serial, atom type and tree identity."""

    def tree(lines):
        return [
            l
            for l in lines
            if l.startswith(("BRANCH", "ENDBRANCH", "TORSDOF", "REMARK INDEX MAP"))
        ]

    source = other.splitlines()
    original = template.splitlines()
    if tree(original) != tree(source):
        raise ValueError("Ligand tree or Meeko serial mapping changed")
    coordinates = {l[6:11]: l for l in source if l.startswith(("ATOM  ", "HETATM"))}
    lines = []
    for l in original:
        if l.startswith(("ATOM  ", "HETATM")):
            update = coordinates.pop(l[6:11])
            if l[77:].strip() != update[77:].strip():
                raise ValueError("Ligand atom type changed")
            l = l[:30] + update[30:54] + l[54:]
        lines.append(l)
    if coordinates:
        raise ValueError("Unmatched ligand atom")
    path.write_text("\n".join(lines) + "\n")


def evaluate(binary, work, output_dir, pdb, meta, file):
    template = (work / pdb / "ligand.pdbqt").resolve()
    receptor = work / pdb / "receptor.pdbqt"
    args = [str(binary), "--receptor", str(receptor), "--seed", "63419", "--cpu", "1"]
    for ax, val in zip("xyz", meta["center"]):
        args += ["--center_" + ax, str(val), "--size_" + ax, "24"]
    relaxed = output_dir / (file.stem + "-relaxed.pdbqt")
    cmd = args + ["--local_only", "--ligand", str(file), "--out", str(relaxed)]
    subprocess.run(cmd, capture_output=True, text=True, check=True)
    rows = []
    for label, pose in (("initial", file), ("relaxed", relaxed)):
        cmd = args + [
            "--score_only",
            "--ligand",
            str(template),
            "--pose_coordinates",
            str(pose),
            "--unbound_energy",
            "0",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)

        def value(key):
            return float(re.search(key + r"\s*:\s*([-+.\d]+)", proc.stdout)[1])

        inter = value("Final Intermolecular Energy")
        intra = value("Final Total Internal Energy")
        rows.append(
            dict(
                pdb=pdb,
                target=meta["target"],
                start=file.stem,
                stage=label,
                rmsd=rms(read_pose(pose, meta["serial"]), np.asarray(meta["mapped"])),
                inter=inter,
                intra=intra,
                total=round(inter + intra, 3),
            )
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("rotation_work", type=Path, help="output of rotation/run.py")
    parser.add_argument("binary", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    rows, failures = [], []
    for folder in sorted((args.rotation_work / "diagonal").iterdir()):
        pdb = folder.name
        case = args.inputs / pdb
        meta = json.loads((case / "metadata.json").read_text())
        template = (case / "ligand.pdbqt").read_text()
        out = args.out / pdb
        out.mkdir(parents=True, exist_ok=True)
        starts = {}
        try:
            starts["projected"] = project(
                template, meta["serial"], np.asarray(meta["mapped"][0])
            )
        except Exception as exc:
            failures.append(dict(pdb=pdb, step="projection", error=str(exc)))
        for seed in (63419, 63421):
            pool = (folder / f"{seed}-geometry.pdbqt").read_text()
            starts[str(seed)] = pool.split("ENDMDL", 1)[0]
        for name, text in starts.items():
            try:
                pose = out / (name + ".pdbqt")
                write_coordinates(template, text, pose)
                rows += evaluate(
                    args.binary.resolve(), args.inputs, out, pdb, meta, pose
                )
            except (subprocess.CalledProcessError, ValueError) as exc:
                failures.append(
                    dict(
                        pdb=pdb,
                        step=name,
                        error=str(getattr(exc, "stderr", None) or exc),
                    )
                )
    with (args.out / "energies.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.out / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")


if __name__ == "__main__":
    main()
