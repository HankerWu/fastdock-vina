"""Post-hoc ring exclusion, rigid closure-distance and production score replay."""

import argparse
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from run import nodes
from analyze import summarize


def main(args):
    data = pd.concat(
        [pd.read_csv(args.archived / "cells.csv"), pd.read_csv(args.s6k / "cells.csv")]
    )
    macro = set()
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        pool = Path(tmp) / "top.pdbqt"
        grid = Path(tmp) / "nodes"
        nodes(grid, 1024)
        samples = np.loadtxt(grid)
        for r in data.query("nodes==1024").itertuples():
            if r.target.startswith("Gal"):
                directory = args.archived / "galectin_affinity"
            elif r.target == "S6K1":
                directory = args.s6k
            else:
                directory = args.archived / "torsion_kinase_affinity"
            ligand = directory / "prepared" / (str(r.ligand) + ".pdbqt")
            if not any(
                re.match(r"^(ATOM|HETATM).*\sG[0-3]\s*$", line)
                for line in ligand.read_text().splitlines()
            ):
                continue
            macro.add((r.target, r.ligand))
            assert r.hinge_moment == 0 and r.joint == r.rigid and abs(r.hinge) < 1e-9
            arm = "pro" if r.target == "S6K1" else "native"
            raw = (
                (directory / "pools" / f"{r.state}_{r.ligand}_{r.seed}_{arm}.pdbqt")
                .read_text()
                .split("ENDMDL", 1)[0]
            )
            pool.write_text(raw + "ENDMDL\n")
            atoms = [
                line for line in raw.splitlines() if line.startswith(("ATOM", "HETATM"))
            ]
            xyz = np.array([[float(s[j : j + 8]) for j in [30, 38, 46]] for s in atoms])
            masks = [
                s.split()[-1] not in ("HD", "H", "G0", "G1", "G2", "G3") for s in atoms
            ]
            center = xyz[masks].mean(0)
            original = pdist(xyz)
            err = max(
                np.max(
                    abs(
                        pdist((xyz - center) @ u[3:12].reshape(3, 3).T + center + u[:3])
                        - original
                    )
                )
                for u in samples
            )
            unbound = float(re.search(r"^REMARK UNBOUND:\s*([-0-9.]+)", raw, re.M)[1])
            cmd = [
                args.pro,
                "--score_only",
                "--ligand",
                ligand,
                "--receptor",
                directory / "prepared" / (r.state + ".pdbqt"),
                "--pose_coordinates",
                pool,
                "--unbound_energy",
                str(unbound),
                "--cpu",
                "1",
            ]
            for k, v in zip("xyz", center):
                cmd += ["--center_" + k, str(v), "--size_" + k, "60"]
            text = subprocess.check_output(cmd, text=True)
            replay = float(
                re.search(r"Estimated Free Energy of Binding\s*:\s*([-0-9.]+)", text)[1]
            )
            rows.append(
                dict(
                    target=r.target,
                    ligand=r.ligand,
                    seed=r.seed,
                    state=r.state,
                    replay=replay,
                    printed=r.native,
                    n_heavy=int(sum(masks)),
                    max_distance_change=err,
                    delta=replay - r.native,
                )
            )
    args.out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out / "macrocycles.csv", index=False)
    keep = [(t, l) not in macro for t, l in zip(data.target, data.ligand)]
    summarize(
        data.loc[keep].query("nodes==1024"),
        pd.read_csv(args.s6k / "docks.csv"),
        args.out / "without_macrocycles",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("archived", "s6k", "pro", "out"):
        p.add_argument(key, type=Path)
    main(p.parse_args())
