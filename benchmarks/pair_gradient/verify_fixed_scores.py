"""Check the full printed Vina energy decomposition on identical archived coordinates."""

import argparse
import csv
import json
import subprocess
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main(binaries, output):
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for source, dest in [
            ("pose/inputs.tar.xz", "inputs"),
            ("basins/starts_and_relaxations.tar.xz", "starts"),
        ]:
            with tarfile.open(HERE.parent / "geometry_comparison" / source) as ar:
                ar.extractall(work / dest, filter="data")
        for case in sorted((work / "starts").iterdir()):
            inputs = work / "inputs" / case.name
            meta = json.loads((inputs / "metadata.json").read_text())
            common = [
                "--score_only",
                "--receptor",
                inputs / "receptor.pdbqt",
                "--ligand",
                inputs / "ligand.pdbqt",
                "--unbound_energy",
                "0",
                "--seed",
                "63419",
                "--cpu",
                "1",
            ]
            for axis, val in zip("xyz", meta["center"]):
                common += ["--center_" + axis, str(val), "--size_" + axis, "24"]
            for pose in sorted(case.glob("*.pdbqt")):
                decompositions = []
                for binary in binaries:
                    log = subprocess.check_output(
                        [binary, *common, "--pose_coordinates", pose], text=True
                    )
                    components = [
                        s
                        for s in log.splitlines()
                        if "(kcal/mol)" in s or "Ligand - " in s or "Flex   - " in s
                    ]
                    if len(components) != 10:
                        raise ValueError(
                            f"Incomplete fixed-coordinate decomposition in {pose}: {components}"
                        )
                    decompositions.append(components)
                if not all(v == decompositions[0] for v in decompositions):
                    raise ValueError(f"Single-point energy mismatch in {pose}")
                rows.append(
                    {
                        "pdb": case.name,
                        "input": pose.stem,
                        "equal_components": len(components),
                        "arms": len(binaries),
                    }
                )
    with output.open("w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binaries", nargs=4, type=Path, help="parent linear cubic shifted")
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    main(a.binaries, a.output)
