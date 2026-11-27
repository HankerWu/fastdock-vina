"""Factorial rigid placement of historical locked-ring outputs, for offline diagnosis only."""

import argparse
import json
import re
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np

EXTERNAL = Path(__file__).resolve().parents[1] / "macrocycle_proposal" / "external"
sys.path.insert(0, str(EXTERNAL))
from external_run import write
from external_score import poses, reference


def orientation(x, y):
    u, _, vt = np.linalg.svd(x.T @ y)
    d = np.eye(3)
    d[-1, -1] = np.linalg.det(u @ vt)
    return u @ d @ vt


def main(a):
    out = []
    with (
        tarfile.open(EXTERNAL / "outputs.tar.xz") as archive,
        tempfile.TemporaryDirectory() as temp,
    ):
        dest = Path(temp) / "pose.pdbqt"
        for pdb in ("5T9W", "7MBO"):
            directory = a.inputs / pdb
            serial, refs, _, _, _, _, tree, _ = reference(a.inputs, pdb)
            center = np.array(
                json.loads((directory / "preparation.json").read_text())["center"]
            )
            box = [
                v
                for i, c in zip("xyz", center)
                for v in ("--center_" + i, str(c), "--size_" + i, "26")
            ]
            for seed in (458921, 740317, 191009):
                name = f"persistent/raw/{pdb}_{seed}_12000_lock.pdbqt"
                lines = (
                    archive.extractfile(name)
                    .read()
                    .decode()
                    .split("ENDMDL")[0]
                    .splitlines()
                )
                dest.write_text("\n".join(lines) + "\n")
                pose = next(poses(dest))
                xyz = np.array([pose[i][1] for i in serial])
                cx = xyz.mean(0)
                x = xyz - cx
                cr = refs[0].mean(0)
                fits = [orientation(x, y - cr) for y in refs]
                best = np.argmin(
                    [np.linalg.norm(x @ r - y + cr) for r, y in zip(fits, refs)]
                )
                rot = fits[best]
                _, axes = np.linalg.eigh(x.T @ x)
                transforms = [
                    ("source", np.eye(3), cx),
                    ("crystal_centroid", np.eye(3), cr),
                    ("crystal_orientation", rot, cx),
                    ("both_crystal", rot, cr),
                ]
                transforms.extend(
                    (f"axis_{i}", 2 * np.outer(v, v) - np.eye(3), cx)
                    for i, v in enumerate(axes.T)
                )
                for label, r, c in transforms:
                    xyznew = x @ r + c
                    rms = float(
                        np.min(np.sqrt(np.square(refs - xyznew).sum(-1).mean(-1)))
                    )
                    formatted = []
                    for line in lines:
                        if line.startswith(("ATOM  ", "HETATM")):
                            v = (pose[int(line[6:11])][1] - cx) @ r + c
                            line = (
                                line[:30] + "".join(f"{z:8.3f}" for z in v) + line[54:]
                            )
                        formatted.append(line)
                    dest.write_text("\n".join(formatted) + "\nENDMDL\n")
                    proc = subprocess.run(
                        [
                            a.parent,
                            "--receptor",
                            directory / "receptor.pdbqt",
                            "--ligand",
                            directory / "ligand.pdbqt",
                            *box,
                            "--pose_coordinates",
                            dest,
                            "--score_only",
                            "--cpu",
                            "1",
                            "--seed",
                            str(seed),
                            "--unbound_energy",
                            "0",
                        ],
                        capture_output=True,
                        check=False,
                        text=True,
                    )
                    values = re.findall(r":\s*([-+.\deE]+)\s+\(kcal/mol\)", proc.stdout)
                    if (
                        len(values) != 10
                        and "outside" not in (proc.stderr + proc.stdout).lower()
                    ):
                        raise ValueError(proc.stdout + proc.stderr)
                    out.append(
                        {
                            "pdb": pdb,
                            "seed": seed,
                            "intervention": label,
                            "rmsd_A": rms,
                            "clash_atoms_2A": sum(tree.query(xyznew)[0] < 2),
                            "centroid_error_A": np.linalg.norm(c - cr),
                            "fixed_zero_ref_score": values[0] if values else "outside",
                        }
                    )
    write(a.output, out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for s in ("inputs", "parent", "output"):
        p.add_argument("--" + s, type=Path, required=True)
    main(p.parse_args())
