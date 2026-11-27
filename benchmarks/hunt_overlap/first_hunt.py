"""Measure the same first Monte Carlo proposal on two energy surfaces."""
import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "pair_gradient/pose"))
import paired


def inspect_pose(prefix, stage, serial, references, receptor, log):
    matches = re.findall(r"GridProbe " + stage + r" (\d+) (\S+) (\S+) (\S+) (\S+)", log)
    grid = np.array(matches, float)
    if len(grid) != len(serial):
        raise ValueError("Inconsistent grid probes and heavy atoms: " + prefix.name)
    raw = grid[:, 1]
    field_grad = np.linalg.norm(grid[:, 2:], axis=1)
    target = paired.atoms(Path(str(prefix) + "." + stage + ".pdbqt").read_text().splitlines())
    xyz = np.array([target[i][0] for i in serial])
    d = cdist(xyz, [item[0] for item in receptor.values()])
    radii = np.array([target[i][1] for i in serial])[:, None] + [item[1] for item in receptor.values()]
    diff = xyz - references
    native = (10 / (10 + np.maximum(raw, 0))) ** 2
    return dict(rmsd_A=float(np.sqrt(np.sum(diff * diff, axis=2).mean(axis=1)).min()),
                clashes=int((d < 0.75 * radii).sum()), min_distance_A=d.min(),
                grid_atoms=len(raw), above_cap=int((raw > 10).sum()), peak_raw_kcal=raw.max(),
                max_tangent_to_native_force_ratio=float(np.max(np.where((raw > 10) & (field_grad > 0),
                                                                         0.25 / native, 1))))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("inputs", type=Path)
    ap.add_argument("probe_binary", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    work = args.inputs.resolve()
    dest = args.output.resolve()
    raw = dest / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    rows = []
    for query in ("1K7L", "2P54", "4UYE", "5FG4", "1A28", "1SQN"):
        receiver = "1A28" if query == "1SQN" else query
        a, b = (json.loads((work / i / "preparation.json").read_text()) for i in (query, receiver))
        rot, shift, _, _ = paired.base.alignment(a, b)
        refs = np.asarray(a["mapped"]) @ rot + shift
        receptor = paired.atoms((work / receiver / "receptor.pdbqt").read_text().splitlines())
        for label, seed in (("random", 89173), ("random", 89179), ("projected", 89173)):
            ligand = (work / query / "ligand.pdbqt" if label == "random" else work / "crystal_diagnostic"
                      / f"{query}-{receiver}-best_projected-input.pdbqt")
            original = None
            for arm in ("native", "all"):
                prefix = raw / f"{query}-{receiver}-{seed}-{label}-{arm}"
                env = {**os.environ, "FASTDOCK_PROBE_OUT": str(prefix)}
                for key in ("FASTDOCK_HUNT_TANGENT", "FASTDOCK_PROBE_PREPARED"):
                    env.pop(key, None)
                if arm == "all":
                    env["FASTDOCK_HUNT_TANGENT"] = arm
                if label == "projected":
                    env["FASTDOCK_PROBE_PREPARED"] = "1"
                cmd = [args.probe_binary.resolve(), "--receptor", work / receiver / "receptor.pdbqt",
                       "--ligand", ligand, "--cpu", "1", "--seed", str(seed), "--exhaustiveness", "1",
                       "--num_modes", "1", "--max_evals", "1", "--search_counters",
                       "--out", str(prefix) + ".output.pdbqt"]
                for axis, center in zip("xyz", b["center"]):
                    cmd.extend(["--center_" + axis, str(center), "--size_" + axis, "24"])
                log = subprocess.run(cmd, env=env, capture_output=True, check=True, text=True).stdout
                Path(str(prefix) + ".log").write_text(log)
                counts = re.search(r"hunt_calls=1 hunt_evals=(\d+) second_calls=1 second_evals=(\d+)", log)
                if not counts:
                    raise ValueError("Expected one hunt and one save")
                for stage in ("start", "end"):
                    measured = inspect_pose(prefix, stage, a["serial"], refs, receptor, log)
                    if stage == "start":
                        coordinates = Path(str(prefix) + ".start.pdbqt").read_bytes()
                        if original and original != (coordinates, measured):
                            raise ValueError("Arms did not start at identical coordinates and field")
                        original = (coordinates, measured)
                    rows.append(dict(query=query, receiver=receiver, start=label, seed=seed, arm=arm, stage=stage,
                                     hunt_evals=int(counts[1]), save_evals=int(counts[2]), **measured))
    paired.write(dest / "rows.csv", rows)


if __name__ == "__main__":
    main()
