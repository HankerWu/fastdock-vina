"""Run matched Late searches sequentially per job; keep raw modes and logs."""

import argparse
import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from freeze import RECEPTORS, SEEDS


def run(work, binary, jobs, arm, seeds=SEEDS):
    with (work / "panel.csv").open() as f:
        panel = list(csv.DictReader(f))
    dest = work / ("pools" if arm == "native" else f"{arm}-pools")
    flags = {
        "native": [],
        "mixed": ["--mixed_vinardo_search"],
        "soft": ["--soft_receptor_search"],
    }[arm]
    dest.mkdir(exist_ok=True)

    def dock_one(case):
        row, pdb, seed = case
        base = dest / f"{pdb}_{row['cid']}_{seed}"
        out = base.with_suffix(".pdbqt")
        if out.exists() and out.read_text().endswith("ENDMDL\n"):
            return
        log = base.with_suffix(".log")
        with log.open("w") as fd:
            subprocess.run(
                [
                    binary,
                    "--receptor",
                    work / "prepared" / f"{pdb}.pdbqt",
                    "--ligand",
                    work / "prepared" / f"{row['cid']}.pdbqt",
                    "--config",
                    work / "prepared" / f"{pdb}.box",
                    "--seed",
                    str(seed),
                    "--cpu",
                    "1",
                    "--exhaustiveness",
                    "2",
                    "--max_evals",
                    "25000",
                    "--num_modes",
                    "9",
                    "--pose_selection",
                    "late",
                    "--energy_range",
                    "100",
                    *flags,
                    "--out",
                    out,
                ],
                stdout=fd,
                stderr=subprocess.STDOUT,
                check=True,
            )

    cases = [
        (row, pdb, seed)
        for seed in seeds
        for row in panel
        for pdb, _ in RECEPTORS[row["target"]]
    ]
    with ThreadPoolExecutor(jobs) as pool:
        for _ in pool.map(dock_one, cases):
            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--jobs", default=3, type=int)
    parser.add_argument("--arm", choices=("native", "mixed", "soft"), default="native")
    parser.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    args = parser.parse_args()
    run(args.work, args.binary.resolve(), args.jobs, args.arm, args.seeds)
