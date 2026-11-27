"""Cross the two AT9283 proposal histories with both native hunt objectives."""

import argparse
import hashlib
import json
import os
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "basin_trace"))
import replay_basin
from summarize import save

GRAPHS = ("inferred", "complete")


def main(work, trace, binary, out):
    out.mkdir(parents=True, exist_ok=False)
    center = json.loads((work / "prepared/2X6E.json").read_text())["center"]
    box = [
        arg
        for axis, value in zip("xyz", center, strict=True)
        for arg in (f"--center_{axis}", str(value), f"--size_{axis}", "26")
    ]
    ligand = work / "prepared/at9283.pdbqt"
    receptor = work / replay_basin.replay.PROTEINS["1OL7_released"]
    bonds = work / "raw/at9283_chemical.bonds"
    costs, cells = [], []

    def run(destination, name, environment):
        stem = out / name
        args = [
            binary,
            "--ligand",
            ligand,
            "--receptor",
            receptor,
            *box,
            "--seed",
            "61231",
            "--exhaustiveness",
            "2",
            "--cpu",
            "1",
            "--num_modes",
            "9",
            "--energy_range",
            "100",
            "--max_evals",
            "25000",
            "--out",
            stem.with_suffix(".pdbqt"),
            "--pose_selection",
            "late",
            "--search_counters",
        ]
        if destination == "complete":
            args += ["--ligand_bonds", bonds]
        env = {
            k: v
            for k, v in os.environ.items()
            if k not in ("FD_DUMP_PROPOSALS", "FD_QUENCH_INPUT", "FD_QUENCH_OUTPUT")
        }
        tic = time.monotonic()
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        with stem.with_suffix(".log").open("w") as f:
            subprocess.run(
                args,
                env=env | {k: str(v) for k, v in environment.items()},
                stdout=f,
                stderr=subprocess.STDOUT,
                check=True,
            )
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        pose = stem.with_suffix(".pdbqt").read_bytes()
        assert pose == (trace / f"at9283_{destination}_61231.pdbqt").read_bytes()
        log = stem.with_suffix(".log").read_text()
        costs.append(
            dict(
                process=name,
                wall_s=time.monotonic() - tic,
                cpu_s=after.ru_stime
                + after.ru_utime
                - before.ru_stime
                - before.ru_utime,
                blind_evals=sum(
                    map(int, re.findall(r"\b(?:hunt|second|explicit)_evals=(\d+)", log))
                ),
                output_sha256=hashlib.sha256(pose).hexdigest(),
            )
        )

    for source in GRAPHS:
        run(source, "capture_" + source, {"FD_DUMP_PROPOSALS": out / (source + ".q")})
    outputs = {}
    for source in GRAPHS:
        for dest in GRAPHS:
            name = source + "_on_" + dest
            values = out / (name + ".xyz")
            run(
                dest,
                name,
                {"FD_QUENCH_INPUT": out / (source + ".q"), "FD_QUENCH_OUTPUT": values},
            )
            outputs[source, dest] = values
    # Only now read experimental ligand coordinates for retrospective assessment.
    _, _, target_order, refs = replay_basin.replay.parameters(work, "at9283")
    for source in GRAPHS:
        obs_order, raw = replay_basin.read_trace(trace / f"at9283_{source}_61231.bin")
        starts = raw[(raw["phase"] == 1) & (raw["local"] == 1) & (raw["event"] == 2)]
        returns = raw[raw["event"] == 4]
        q = np.loadtxt(out / (source + ".q"))
        assert q.shape == (len(starts), 14) and len(starts) == len(returns)
        assert np.array_equal(
            q[:, :2], np.stack([starts["traj"], starts["step"]], axis=1)
        )
        assert np.array_equal(
            q[:, :2], np.stack([returns["traj"], returns["step"]], axis=1)
        )
        d = {}
        for dest in GRAPHS:
            values = outputs[source, dest]
            with values.open() as f:
                serials = list(map(int, f.readline().removeprefix("#").split()))
            assert serials == obs_order
            results = np.loadtxt(values)
            assert np.array_equal(results[:, :2], q[:, :2])
            xyz = results[:, 4:].reshape(len(q), len(serials), 3)
            reindexed = refs[:, [target_order.index(s) for s in serials]]

            def rmsd(x):
                return np.sqrt(np.square(x[:, None] - reindexed).sum(-1).mean(-1)).min(
                    axis=1
                )

            if dest == source:
                assert np.array_equal(xyz, returns["xyz"])
                assert np.array_equal(results[:, 3], returns["energy"])
            entry = dict(
                serials=serials,
                generalized=q[:, 2:],
                trajectory=q[:, 0].astype(int),
                step=q[:, 1].astype(int),
                start_xyz=starts["xyz"],
                return_xyz=xyz,
                evals=results[:, 2].astype(int),
                energy=results[:, 3],
                start_rmsd=rmsd(starts["xyz"]),
                end_rmsd=rmsd(xyz),
            )
            np.savez_compressed(out / (source + "_on_" + dest + ".npz"), **entry)
            d[dest] = entry
        for i in range(len(q)):
            cells.append(
                dict(
                    source=source,
                    trajectory=int(q[i, 0]),
                    step=int(q[i, 1]),
                    start_A=d[source]["start_rmsd"][i],
                    **{f"end_{g}_A": d[g]["end_rmsd"][i] for g in GRAPHS},
                    **{f"evals_{g}": int(d[g]["evals"][i]) for g in GRAPHS},
                )
            )
    save(out / "paired.csv", cells)
    save(out / "costs.csv", costs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("work", "trace", "binary", "out"):
        p.add_argument(arg, type=Path)
    main(*(v.resolve() for v in vars(p.parse_args()).values()))
