"""Verify ranking windows, a second real receptor and explicitly requested radii."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from compare import models, output


def retained(pool, radius):
    selected = []
    for n, (_, x) in enumerate(pool):
        if all(
            np.sqrt(np.square(x - pool[j][1]).sum(axis=1).mean()) >= radius
            for j in selected
        ):
            selected.append(n)
    return selected


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("gal", "archive4", "ar", "out", "old", "new"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    donor = (
        next(a.archive4.glob("Gal3_3T1L_26_391781_candidate_3000.pdbqt"))
        .read_text()
        .split("ENDMDL")[0]
    )
    elements = {"NA": "N", "OA": "O", "A": "C"}
    coords = []
    for line in donor.splitlines():
        if line.startswith(("HETATM", "ATOM  ")) and line.split()[-1] not in (
            "HD",
            "H",
        ):
            element = line.split()[-1]
            coords.append(
                " ".join(
                    [
                        elements.get(element, element),
                        *[line[i : i + 8].strip() for i in (30, 38, 46)],
                    ]
                )
            )
    template = a.out / "reference.txt"
    template.write_text("\n".join(coords) + "\n")
    gallery = [
        "--config",
        a.gal / "3T1L.box",
        "--receptor",
        a.gal / "3T1L.pdbqt",
        "--ligand",
        a.gal / "10.pdbqt",
        "--seed",
        "391781",
        "--max_evals",
        "3000",
    ]
    ligand = a.ar / "1T7T"
    rec = [
        "--receptor",
        ligand / "receptor.pdbqt",
        "--ligand",
        ligand / "ligand.pdbqt",
        "--seed",
        "438013",
        "--max_evals",
        "12000",
    ]
    for axis, v in zip(
        "xyz", json.loads((ligand / "preparation.json").read_text())["center"]
    ):
        rec.extend(["--center_" + axis, str(v), "--size_" + axis, "22.5"])
    results = []
    cases = [
        (
            "custom" + str(w),
            gallery + ["--pose_template", template, "--template_stage", "contact_rank"],
            1,
            w,
        )
        for w in [100, 3, 1, 0.25, 0]
    ]
    cases += [("native-r" + str(r), gallery, r, 100) for r in [0, 2]] + [
        ("1T7T", rec, 1, 100)
    ]
    for name, common, radius, window in cases:
        last = None
        for arm in ("old", "new"):
            root = a.out / (name + "_" + arm)
            pdb = a.out / (root.name + ".pdbqt")
            log = a.out / (root.name + ".log")
            cmd = [
                getattr(a, arm),
                *common,
                "--out",
                pdb,
                "--cpu",
                "1",
                "--exhaustiveness",
                "4",
                "--num_modes",
                "9",
                "--energy_range",
                str(window),
                "--min_rmsd",
                str(radius),
                "--search_counters",
            ]
            with log.open("w") as f:
                subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT)
            pool = models(pdb.read_bytes())
            counts = "\n".join(
                x
                for x in log.read_text().splitlines()
                if x.startswith(
                    ("Search trajectory:", "Search counters:", "Hunt chain:")
                )
            )
            if arm == "old":
                last = pool, counts
                continue
            old, c0 = last
            idx = retained(old, radius)
            if [x[0] for x in pool] != [old[i][0] for i in idx] or counts != c0:
                raise ValueError(name + " geometry/window mismatch")
            results.append(
                dict(
                    case=name,
                    radius=radius,
                    energy_window=window,
                    old=len(old),
                    new=len(pool),
                    original_ranks=",".join(str(i + 1) for i in idx),
                    counter_sha256=hashlib.sha256(counts.encode()).hexdigest(),
                )
            )
    output(a.out / "edge_cases.csv", results)
