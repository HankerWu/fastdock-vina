"""Replay original, placed, and once-minimized geometries at matched RNG draws."""
import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist, pdist

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "pair_gradient/pose"))
import paired


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "probe_binary", "out"):
        p.add_argument(name, type=Path)
    args = p.parse_args()
    raw = args.out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    rows = []
    for query in ("1K7L", "2P54", "4UYE", "5FG4", "1A28", "1SQN"):
        receiver = "1A28" if query == "1SQN" else query
        a, b = (json.loads((args.inputs / k / "preparation.json").read_text()) for k in (query, receiver))
        rot, shift, _, _ = paired.base.alignment(a, b)
        refs = np.asarray(a["mapped"]) @ rot + shift
        receptor_file = args.inputs / receiver / "receptor.pdbqt"
        rec = paired.atoms(receptor_file.read_text().splitlines())
        for seed in (89173, 89179):
            baseline = None
            for arm in ("native", "fit", "invalid"):
                name = raw / f"{query}-{receiver}-{seed}-{arm}"
                env = {**os.environ, "FASTDOCK_PROBE_OUT": str(name)}
                env.pop("FASTDOCK_FIRST_PLACEMENT", None)
                if arm != "native":
                    env["FASTDOCK_FIRST_PLACEMENT"] = arm
                cmd = [args.probe_binary, "--receptor", receptor_file, "--ligand", args.inputs / query / "ligand.pdbqt",
                       "--cpu", "1", "--seed", str(seed), "--exhaustiveness", "1", "--num_modes", "1",
                       "--max_evals", "1", "--search_counters", "--out", str(name) + ".output.pdbqt"]
                for axis, x in zip("xyz", b["center"]):
                    cmd.extend(["--center_" + axis, str(x), "--size_" + axis, "24"])
                log = subprocess.run(cmd, env=env, capture_output=True, text=True, check=True).stdout
                Path(str(name) + ".log").write_text(log)
                counts = re.search(r"hunt_calls=1 hunt_evals=(\d+) second_calls=1 second_evals=(\d+)", log)
                if not counts:
                    raise ValueError("Expected one hunt and one second local call")
                for stage in ("original", "placed", "end"):
                    atoms = paired.atoms(Path(str(name) + "." + stage + ".pdbqt").read_text().splitlines())
                    xyz = np.array([atoms[i][0] for i in a["serial"]])
                    if stage == "original":
                        if baseline is None:
                            baseline = xyz
                        np.testing.assert_array_equal(xyz, baseline)
                    distortion = np.max(np.abs(pdist(xyz) - pdist(baseline)))
                    if stage == "placed" and distortion > 0.0035:
                        raise ValueError("Placement distorted the input geometry")
                    # Coordinates are emitted to 0.001 Å; allow for this rounding.
                    excess = np.maximum(np.abs(xyz - b["center"]) - 12, 0)
                    if stage == "placed" and arm != "native" and excess.max() > 0.00051:
                        raise ValueError("Unexpected exterior atom after mapping in these fixtures")
                    delta = xyz - refs
                    d = cdist(xyz, [v[0] for v in rec.values()])
                    radii = np.array([atoms[i][1] for i in a["serial"]])[:, None] + [v[1] for v in rec.values()]
                    rows.append(dict(query=query, receiver=receiver, seed=seed, arm=arm, stage=stage,
                                     exterior_atoms=int((excess.max(1) > 0.00051).sum()),
                                     maximum_excess_A=excess.max(), clashes=int((d < 0.75*radii).sum()),
                                     rmsd_A=float(np.sqrt(np.sum(delta*delta,axis=2).mean(axis=1)).min()),
                                     largest_internal_distance_change_A=distortion, hunt_evals=int(counts[1]),
                                     save_evals=int(counts[2]), xyz_unchanged_vs_original=np.array_equal(xyz,baseline)))
    paired.write(args.out / "rows.csv", rows)


if __name__ == "__main__":
    main()
