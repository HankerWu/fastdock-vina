"""Run the prereported thrombin settings with half-strength receptor repulsion in even trajectories."""

import argparse
import itertools
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

RECEIVERS = ("3D49", "2ZFF", "2ZC9")
COMPOUNDS = ("4", "6a", "6e", "6b", "1c", "1b", "3a", "1d", "7a", "1a", "3b", "5")
SEEDS = (83519, 194027)


def main(binary, baseline, work, jobs):
    work.mkdir(parents=True, exist_ok=True)
    for sub in ["raw", "prepared"]:
        if not (work / sub).exists():
            (work / sub).symlink_to(
                (baseline / sub).resolve(), target_is_directory=True
            )
    shutil.copyfile(baseline / "series.csv", work / "series.csv")
    (work / "results").mkdir(exist_ok=True)

    def dock(item):
        rec, i, s = item
        dest = work / "results" / f"{rec}_{i}_{s}.pdbqt"
        with dest.with_suffix(".log").open("w") as h:
            subprocess.run(
                [
                    str(binary),
                    "--config",
                    str(work / "prepared" / f"{rec}.box"),
                    "--receptor",
                    str(work / "prepared" / f"{rec}.pdbqt"),
                    "--ligand",
                    str(work / "prepared" / f"{i}.pdbqt"),
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
                    "--soft_receptor_search",
                    "--out",
                    str(dest),
                ],
                stdout=h,
                stderr=subprocess.STDOUT,
                check=True,
            )

    with ThreadPoolExecutor(jobs) as executor:
        list(executor.map(dock, itertools.product(RECEIVERS, COMPOUNDS, SEEDS)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("--jobs", type=int, default=3)
    a = parser.parse_args()
    main(a.binary.resolve(), a.baseline, a.work, a.jobs)
