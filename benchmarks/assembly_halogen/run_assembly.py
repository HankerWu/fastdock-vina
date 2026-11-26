"""Rebuild biological units and run matched ASU/assembly searches."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

from build_receptor_assembly import build

PROBES = ("7LZQ_YJV", "8DW5_FQ7")
CONTROLS = ("7BJ6_TVK", "7RWO_7WN")
INITIAL = (314159, 271828, 161803)
REPEAT = (142857, 173205, 223607, 244949, 264575, 282843)


def run(work, out, binary, *, high=False):
    jobs = []
    for case in (PROBES + ("2ROX_T44",)) if high else PROBES + CONTROLS:
        folder = work / case
        dest = out / case
        dest.mkdir(exist_ok=True, parents=True)
        meta = build(
            folder / "deposition.pdb",
            folder / "receptor.pdbqt",
            "1",
            dest / "receptor.pdbqt",
        )
        (dest / "assembly.json").write_text(json.dumps(meta, indent=2) + "\n")
        for seed in (
            (316227, 331663, 346410)
            if high
            else INITIAL + (REPEAT if case in PROBES else ())
        ):
            for arm in ("asu", "assembly"):
                if not high and seed in INITIAL and arm == "asu":
                    continue
                receptor = (folder if arm == "asu" else dest) / "receptor.pdbqt"
                for cap in (60000, 240000) if high else (60000,):
                    output = dest / (
                        f"{arm}-{seed}-cap{cap}.pdbqt"
                        if high
                        else f"{arm}-{seed}.pdbqt"
                    )
                    jobs.append((folder, receptor, output, seed, cap))

    def dock(job):
        folder, receptor, output, seed, cap = job
        if output.exists() and output.read_text().endswith("ENDMDL\n"):
            return
        with output.with_suffix(".log").open("w") as f:
            subprocess.run(
                [
                    binary,
                    "--receptor",
                    receptor,
                    "--ligand",
                    folder / "query.pdbqt",
                    "--config",
                    folder / "box",
                    "--seed",
                    str(seed),
                    "--cpu",
                    "1",
                    "--exhaustiveness",
                    "3",
                    "--max_evals",
                    str(cap),
                    "--num_modes",
                    "9",
                    "--pose_selection",
                    "late",
                    "--refine_candidates",
                    "12",
                    "--energy_range",
                    "100",
                    "--out",
                    output,
                ],
                stdout=f,
                stderr=subprocess.STDOUT,
                check=True,
            )

    with ThreadPoolExecutor(3) as pool:
        list(pool.map(dock, jobs))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--high", action="store_true")
    a = p.parse_args()
    run(a.work, a.out, a.binary.resolve(), high=a.high)
