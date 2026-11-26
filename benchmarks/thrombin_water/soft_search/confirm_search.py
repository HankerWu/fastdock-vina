"""Two new fixed seeds on the published thrombin panel. Run both hard/soft search with the same binary."""

import argparse
import itertools
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from soft_challenge import COMPOUNDS, RECEIVERS

CONFIRMATION_SEEDS = (88091, 302123)


def main(binary, root, dest, jobs):
    dest.mkdir(parents=True, exist_ok=True)

    def run(task):
        r, lig, s, arm = task
        out = dest / f"{r}_{lig}_{s}_{arm}.pdbqt"
        args = [
            str(binary),
            "--config",
            str(root / "prepared" / f"{r}.box"),
            "--receptor",
            str(root / "prepared" / f"{r}.pdbqt"),
            "--ligand",
            str(root / "prepared" / f"{lig}.pdbqt"),
            "--cpu",
            "1",
            "--pose_selection",
            "late",
            "--exhaustiveness",
            "4",
            "--max_evals",
            "50000",
            "--num_modes",
            "12",
            "--refine_candidates",
            "16",
            "--energy_range",
            "8",
            "--seed",
            str(s),
            "--out",
            str(out),
            "--search_counters",
        ]
        if arm == "soft":
            args.append("--soft_receptor_search")
        with out.with_suffix(".log").open("w") as h:
            subprocess.run(args, stdout=h, stderr=subprocess.STDOUT, check=True)

    with ThreadPoolExecutor(jobs) as pool:
        list(
            pool.map(
                run,
                itertools.product(
                    RECEIVERS, COMPOUNDS, CONFIRMATION_SEEDS, ("hard", "soft")
                ),
            )
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for name in ["binary", "baseline", "dest"]:
        p.add_argument(name, type=Path)
    p.add_argument("--jobs", type=int, default=3)
    a = p.parse_args()
    main(a.binary.resolve(), a.baseline.resolve(), a.dest.resolve(), a.jobs)
