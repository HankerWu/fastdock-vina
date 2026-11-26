"""Run a matched single- vs all-assignment global search, preserving the original Vina budget."""

import argparse
import itertools
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from prepare import RECEPTORS


def run_one(binary, prepared, work, cell):
    receptor, ligand, seed, arm = cell
    stem = work / f"{receptor}_{ligand}_{seed}_{arm}"
    command = [
        binary,
        "--ligand",
        prepared / f"{ligand}_ligand.pdbqt",
        "--receptor",
        prepared / f"{receptor}_receptor.pdbqt",
        "--config",
        prepared / f"{receptor}_box.txt",
        "--seed",
        str(seed),
        "--exhaustiveness",
        "4",
        "--cpu",
        "1",
        "--max_evals",
        "50000",
        "--pose_selection",
        "late",
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--search_counters",
        "--out",
        stem.with_suffix(".pdbqt"),
    ]
    if arm != "native":
        command.extend(["--search_anchor", work / f"{receptor}_{ligand}.{arm}"])
    with stem.with_suffix(".log").open("w") as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("panel", choices=RECEPTORS)
    p.add_argument("prepared", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--seeds", type=int, nargs="+", default=[314159, 271828])
    p.add_argument(
        "--arms",
        nargs="+",
        choices=("native", "single", "multi"),
        default=["single", "multi"],
    )
    a = p.parse_args()
    compounds = sorted(json.loads((a.work / "graphs.json").read_text()))
    cells = itertools.product(RECEPTORS[a.panel], compounds, a.seeds, a.arms)
    with ThreadPoolExecutor(a.jobs) as pool:
        list(
            pool.map(
                lambda c: run_one(a.binary.resolve(), a.prepared, a.work, c), cells
            )
        )
