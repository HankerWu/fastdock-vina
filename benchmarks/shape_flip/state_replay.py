"""At the first changed proposal, compare both moves before and after the same grid hunt."""

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

EXTERNAL = Path(__file__).resolve().parents[1] / "macrocycle_proposal" / "external"
sys.path.insert(0, str(EXTERNAL))
from external_run import sha256, write
from external_score import poses, reference


def main(a):
    destination = a.output.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    metrics, events = [], []
    refs = {p: reference(a.inputs, p) for p in ("7MBO", "5T9W")}
    for r in csv.DictReader((a.base / "runs.csv").open()):
        if r["arm"] not in ("always", "adaptive") or r["cap"] != "12000":
            continue
        d = a.inputs / r["pdb"]
        meta = json.loads((d / "preparation.json").read_text())
        box = [
            v
            for i, c in zip("xyz", meta["center"])
            for v in ("--center_" + i, str(c), "--size_" + i, "26")
        ]
        key = r["key"]
        prefix = destination / key
        output = destination / "final.pdbqt"
        subprocess.run(
            [
                getattr(a, r["arm"]),
                "--receptor",
                d / "receptor.pdbqt",
                "--ligand",
                d / "ligand.pdbqt",
                *box,
                "--seed",
                r["seed"],
                "--max_evals",
                r["cap"],
                "--exhaustiveness",
                "1",
                "--cpu",
                "1",
                "--num_modes",
                "5",
                "--verbosity",
                "0",
                "--energy_range",
                "100",
                "--search_counters",
                "--out",
                output,
            ],
            env={**os.environ, "SHAPE_TRACE": str(prefix)},
            capture_output=True,
            check=True,
        )
        if sha256(output) != r["sha256"]:
            raise ValueError("Observer changed output " + key)
        output.unlink()
        if int(r["flips"]) == 0:
            continue
        details = list(
            csv.DictReader(prefix.with_suffix(".tsv").open(), delimiter="\t")
        )
        events.extend(dict(key=key, **event) for event in details)
        one = refs[r["pdb"]]
        serial, mapped, sdf_order, sdf_refs, _, _, tree, _ = one
        stages = list(poses(prefix.with_suffix(".pdbqt")))
        if len(stages) != 5:
            raise ValueError(key)
        for stage, pose in zip(
            ("source", "native_proposal", "native_hunt", "flip_proposal", "flip_hunt"),
            stages,
        ):
            xyz = np.array([pose[i][1] for i in serial])
            xyz2 = np.array([pose[i][1] for i in sdf_order])
            rms = float(np.min(np.sqrt(((mapped - xyz) ** 2).sum(-1).mean(-1))))
            independent = float(
                np.min(np.sqrt(((sdf_refs - xyz2) ** 2).sum(-1).mean(-1)))
            )
            if abs(rms - independent) > 0.0015:
                raise ValueError(key)
            # The CLI uses ceil(26/0.375) voxels along each dimension.
            half_span = np.ceil(26 / 0.375) * 0.375 / 2
            outside = np.maximum(
                np.abs(xyz - np.asarray(meta["center"])) - half_span, 0
            )
            metrics.append(
                {
                    "key": key,
                    "stage": stage,
                    "rmsd_A": rms,
                    "sdf_rmsd_A": independent,
                    "summed_wall_violation_A": outside.sum(),
                    "centroid_error_A": np.linalg.norm(xyz.mean(0) - mapped[0].mean(0)),
                    "heavy_clashes_2A": int(np.sum(tree.query(xyz)[0] < 2)),
                }
            )
    write(destination / "stages.csv", metrics)
    write(destination / "decisions.csv", events)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for s in ("inputs", "base", "output", "always", "adaptive"):
        p.add_argument("--" + s, type=Path, required=True)
    main(p.parse_args())
