"""Check paired budgets, rigid geometry, post-relaxation coordinates and native local parity."""

import argparse
import csv
import json
import subprocess
import tempfile
from pathlib import Path

from scipy.spatial.distance import pdist

from run import parse_pose


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    p = argparse.ArgumentParser()
    for name in ("dev", "hold", "full", "vina"):
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    rows = read(a.full / "cost.csv")
    cells = read(a.full / "cells.csv")
    assert len(cells) == 40
    grouped = {}
    for row in rows:
        if row["arm"] in ("lattice", "random"):
            grouped.setdefault((row["query"], row["receiver"], row["seed"]), {})[
                row["arm"]
            ] = row
    assert len(grouped) == 20
    for g in grouped.values():
        for field in ("probes", "lookups", "refines"):
            assert g["lattice"][field] == g["random"][field]
    devmeta = json.loads((a.dev / "metadata.json").read_text())
    max_distortion = 0
    for (q, r, s), _ in grouped.items():
        dev = r in ("5C7N", "4OD0")
        base = a.dev if dev else a.hold
        query = (
            devmeta[q]
            if dev
            else json.loads((base / q / "preparation.json").read_text())
        )
        original = pdist(parse_pose(base / q / "ligand.pdbqt", query["serial"]))
        for arm in ("lattice", "random"):
            for path in (a.full / "runs" / f"{r}_{q}_{s}_{arm}").glob("*.pdbqt"):
                err = abs(pdist(parse_pose(path, query["serial"])) - original).max()
                max_distortion = max(err, max_distortion)
    assert max_distortion < 0.002
    print(
        "all 20 probe, atom-lookup and refinement count pairs match; max pair-distance drift (A)",
        max_distortion,
    )
    for q, r, s, k in [
        ("5EPS", "5C7N", "43103", 4),
        ("5EPS", "5C7N", "43117", 1),
        ("3WKE", "4OD0", "43117", 1),
        ("3WKE", "4OD0", "43117", 3),
    ]:
        path = a.full / "runs" / f"{r}_{q}_{s}_lattice" / f"rank{k}.pdbqt"
        with tempfile.TemporaryDirectory() as w:
            result = Path(w) / "native.pdbqt"
            cmd = [
                str(a.vina),
                "--local_only",
                "--receptor",
                str(a.dev / r / "receptor.pdbqt"),
                "--ligand",
                str(path),
                "--out",
                str(result),
                "--cpu",
                "1",
            ]
            for axis, val in zip("xyz", devmeta[r]["center"]):
                cmd += ["--center_" + axis, str(val), "--size_" + axis, "24"]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
            assert result.read_bytes() == Path(str(path) + ".refined").read_bytes()
    print(
        "four native --local_only outputs byte-identical to observer, including both competitors in 3WKE."
    )


if __name__ == "__main__":
    main()
