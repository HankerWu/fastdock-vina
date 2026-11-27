"""Paired native AutoDock4-map and Vina fast-force-field search on the FKBP12 ester series."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

RECEPTORS = ("1FKG", "1FKH")
SEEDS = (131071, 196613, 262147)
COMPOUNDS = range(4, 12)
TYPES = ("A", "C", "N", "OA")
SPACING = 0.375


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def center(work, receiver):
    xyz = json.loads((work / receiver / "preparation.json").read_text())["center"]
    # Match the center precision in AutoGrid's ASCII map header.
    return tuple(float(f"{x:.3f}") for x in xyz)


def prep_grid(work, out, autogrid, receiver):
    dest = out / "maps" / receiver
    dest.mkdir(parents=True, exist_ok=True)
    receptor = work / receiver / "receptor.pdbqt"
    loc = dest / "receptor.pdbqt"
    if loc.exists():
        if sha(loc) != sha(receptor):
            raise ValueError("receptor input changed")
    else:
        os.symlink(receptor, loc)
    rtypes = sorted(
        {
            line[77:].strip()
            for line in receptor.read_text().splitlines()
            if line.startswith(("ATOM", "HETATM"))
        }
    )
    grid = [
        "npts 64 64 64",
        "gridfld receptor.maps.fld",
        "spacing 0.375",
        "receptor_types " + " ".join(rtypes),
        "ligand_types " + " ".join(TYPES),
        "receptor receptor.pdbqt",
        "gridcenter " + " ".join(map(str, center(work, receiver))),
        "smooth 0.500",
    ]
    grid.extend("map receptor." + t + ".map" for t in TYPES)
    grid.extend(
        ["elecmap receptor.e.map", "dsolvmap receptor.d.map", "dielectric -42.000"]
    )
    (dest / "receptor.gpf").write_text("\n".join(grid) + "\n")
    subprocess.run(
        [autogrid, "-p", "receptor.gpf", "-l", "receptor.glg"],
        cwd=dest,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    if "Successful Completion" not in (dest / "receptor.glg").read_text():
        raise ValueError("Incomplete AutoGrid map construction")
    return {p.name: sha(p) for p in sorted(dest.glob("*.map"))}


def run_one(work, ligands, out, vina, job, *, search_counters=False):
    receiver, compound, seed, arm = job
    stem = out / "pools" / f"{compound}_{receiver}_{seed}_{arm}"
    pose = stem.with_suffix(".pdbqt")
    log = stem.with_suffix(".log")
    args = [
        str(vina),
        "--ligand",
        str(ligands / f"{compound}.pdbqt"),
        "--scoring",
        arm,
        "--cpu",
        "1",
        "--exhaustiveness",
        "3",
        "--max_evals",
        "24000",
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--seed",
        str(seed),
        "--out",
        str(pose),
    ]
    if search_counters:
        args.append("--search_counters")
    if arm == "ad4":
        args += ["--maps", str(out / "maps" / receiver / "receptor")]
    else:
        args += ["--receptor", str(work / receiver / "receptor.pdbqt")]
        for axis, value in zip("xyz", center(work, receiver)):
            args += ["--center_" + axis, str(value), "--size_" + axis, "24"]
    start = time.monotonic()
    with log.open("w") as f:
        subprocess.run(
            args, stdout=f, stderr=subprocess.STDOUT, check=True, timeout=900
        )
    duration = time.monotonic() - start
    text = pose.read_text()
    scores = [
        float(s.split()[3])
        for s in text.splitlines()
        if s.startswith("REMARK VINA RESULT:")
    ]
    if not text.rstrip().endswith("ENDMDL") or len(scores) != text.count("ENDMDL"):
        raise ValueError(f"Incomplete {pose}")
    return dict(
        receiver=receiver,
        compound=compound,
        seed=seed,
        arm=arm,
        score=scores[0],
        modes=len(scores),
        seconds=duration,
        sha256=sha(pose),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "ligands", "out", "autogrid", "vina"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    a = parser.parse_args()
    work, ligands, out, ag, vina = [
        getattr(a, n).resolve() for n in ("work", "ligands", "out", "autogrid", "vina")
    ]
    (out / "pools").mkdir(parents=True, exist_ok=True)
    for compound in COMPOUNDS:
        types = {
            l[77:].strip()
            for l in (ligands / f"{compound}.pdbqt").read_text().splitlines()
            if l.startswith(("ATOM", "HETATM"))
        }
        if not types.issubset(TYPES):
            raise ValueError(
                f"Missing grid for ligand {compound}: {types - set(TYPES)}"
            )
    with ThreadPoolExecutor(2) as exe:
        maps = dict(
            zip(RECEPTORS, exe.map(lambda r: prep_grid(work, out, ag, r), RECEPTORS))
        )
    metadata = {
        "maps": maps,
        "ligands": {str(n): sha(ligands / f"{n}.pdbqt") for n in COMPOUNDS},
        "centers_A": {r: center(work, r) for r in RECEPTORS},
        "receptors": {r: sha(work / r / "receptor.pdbqt") for r in RECEPTORS},
        "autogrid_sha256": sha(ag),
        "vina_sha256": sha(vina),
        "seeds": SEEDS,
        "per_trajectory_max_evals": 24000,
        "exhaustiveness": 3,
        "cpu": 1,
        "npts": [64] * 3,
        "spacing_A": SPACING,
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    jobs = itertools.product(RECEPTORS, COMPOUNDS, SEEDS, ("vina", "ad4"))
    with ThreadPoolExecutor(a.jobs) as exe:
        rows = list(exe.map(lambda j: run_one(work, ligands, out, vina, j), jobs))
    table(out / "cells.csv", rows)


if __name__ == "__main__":
    main()
