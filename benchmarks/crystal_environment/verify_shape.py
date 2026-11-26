"""Replay both fixed donor conditions on every pose using the production C++ header."""

import argparse
import itertools
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from evaluate import CHAIN, ROOT, models, read
from packing import atoms
from shape import transfer


def format_atoms(xyz, elements):
    return "".join(
        f"{e.capitalize()} {x:.10f} {y:.10f} {z:.10f}\n"
        for (x, y, z), e in zip(xyz, elements, strict=True)
    )


def run(a):
    paths = {
        "hcaii": (a.ca, a.ca_runs),
        "brd4": (a.brd, a.brd_runs),
        "t4l": (a.t4, a.t4_runs),
    }
    expected = []
    inputs = []
    data = read(a.results / "shape_terms.csv")
    for target, (source, runs) in paths.items():
        name = "metadata.json" if target == "hcaii" else "manifest.json"
        meta = json.loads((source / "prepared" / name).read_text())
        rows = [r for r in data if r["target"] == target]
        keys = lambda r: (r["receptor"], r["ligand"], r["seed"], r["slot"], r["donor"])
        for (rec, lig, seed, slot, donor), group in itertools.groupby(
            sorted(rows, key=keys), key=keys
        ):
            assert donor != lig
            poses = models(
                runs / f"{rec}_{lig}_{seed}_{slot}.pdbqt",
                None if target == "hcaii" else meta[lig]["heavy_serials"],
            )
            component = meta[donor].get("ccd") or meta[donor]["component"]
            pdb = source / "raw" / f"{donor}.pdb"
            reference = [
                atom
                for atom in atoms(pdb)
                if atom.residue == component
                and atom.chain == CHAIN.get(donor, "A")
                and (target != "hcaii" or atom.number == "303")
            ]
            origin, orientation, destination = transfer(
                source / "raw" / f"{rec}.pdb", pdb, rec, donor, target
            )
            donor_text = format_atoms(
                [a.xyz for a in reference], [a.element for a in reference]
            )
            for r in group:
                expected.append(float(r["existing"]))
                xyz, elems = poses[int(r["model"])]
                inputs.extend(
                    [
                        f"{len(reference)} {len(elems)}\n",
                        donor_text,
                        format_atoms((xyz - origin) @ orientation + destination, elems),
                    ]
                )
    with tempfile.TemporaryDirectory() as folder:
        binary = Path(folder) / "compare"
        subprocess.run(
            [
                "c++",
                "-std=c++17",
                "-O2",
                "-I",
                str(ROOT / "src/lib"),
                str(Path(__file__).with_name("compare_template.cpp")),
                "-o",
                str(binary),
            ],
            check=True,
        )
        result = subprocess.run(
            [str(binary)],
            input="".join(inputs),
            text=True,
            capture_output=True,
            check=True,
        )
    output = np.fromstring(result.stdout, sep="\n")
    if len(expected) != len(output) or len(output) != len(data):
        raise ValueError("Incomplete output")
    err = np.max(np.abs(expected - output))
    if err > 1e-9:
        raise ValueError(err)
    print(f"{len(output)} pose/donor overlaps matched original C++ to {err:.3g}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("ca", "brd", "t4", "ca_runs", "brd_runs", "t4_runs", "results"):
        p.add_argument(arg, type=Path)
    run(p.parse_args())
