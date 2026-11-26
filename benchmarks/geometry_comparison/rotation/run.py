"""Run the fixed isotropic rotation control on the measured cognate panel."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "geometry_metric"))
import panel


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "inputs", type=Path, help="fresh pose_prepare.py output directory"
    )
    parser.add_argument("work", type=Path)
    parser.add_argument("diagonal", type=Path)
    parser.add_argument("isotropic", type=Path)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    with (HERE.parent / "pose_cases.csv").open() as f:
        cases = [r for r in csv.DictReader(f) if r["cohort"] != "reciprocal_candidate"]
    measures = []
    for row in cases:
        pdb = row["pdb"]
        inputs = (args.inputs / pdb).resolve()
        if not (inputs / "metadata.json").exists():
            continue
        meta = json.loads((inputs / "metadata.json").read_text())
        assert meta["pdb"] == pdb and meta["status"] == "ready"
        for i, seed in enumerate(panel.SEEDS):
            arms = [("diagonal", args.diagonal), ("isotropic", args.isotropic)]
            for name, binary in arms[:: 1 if i == 0 else -1]:
                out = args.work / name / pdb
                out.mkdir(parents=True, exist_ok=True)
                for file in ("ligand.pdbqt", "receptor.pdbqt"):
                    location = out / file
                    if not location.exists():
                        location.symlink_to(inputs / file)
                modes = panel.dock(out.parent, binary.resolve(), meta, "geometry", seed)
                measures.append(
                    dict(
                        pdb=pdb,
                        target=meta["target"],
                        seed=seed,
                        arm=name,
                        top1_rmsd=modes[0]["rmsd"],
                        best_rmsd=min(m["rmsd"] for m in modes),
                        top1_score=modes[0]["score"],
                        modes=len(modes),
                        output_sha256=digest(out / f"{seed}-geometry.pdbqt"),
                    )
                )
        print(pdb, flush=True)
    with (args.work / "results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(measures[0]))
        writer.writeheader()
        writer.writerows(measures)


if __name__ == "__main__":
    main()
