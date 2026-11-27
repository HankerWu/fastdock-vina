"""Replay the fixed local starts and score their coordinates on the parent ligand graph."""

import argparse
import csv
import itertools
import json
import re
import subprocess
import sys
import tarfile
from pathlib import Path
from time import perf_counter

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "paired_crossdocking"))
from crystal_start import read_pose, rms


def call(argv):
    return subprocess.run(argv, check=True, text=True, capture_output=True).stdout


def common_score(binary, args, ligand, coordinates):
    log = call(
        [
            binary,
            *args,
            "--score_only",
            "--ligand",
            ligand,
            "--pose_coordinates",
            coordinates,
            "--unbound_energy",
            "0",
        ]
    )

    def component(label):
        return float(re.search(label + r"\s*:\s*([-+.\d]+)", log)[1])

    return round(
        component("Final Intermolecular Energy")
        + component("Final Total Internal Energy"),
        3,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("parent", "linear", "cubic", "out"):
        parser.add_argument(name, type=Path)
    parser.add_argument(
        "--shift", type=Path, help="Include exploratory potential-shift arm"
    )
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    source = HERE.parent / "geometry_comparison"
    for archive, dest in [
        ("pose/inputs.tar.xz", "fixtures"),
        ("basins/starts_and_relaxations.tar.xz", "starts"),
    ]:
        with tarfile.open(source / archive) as stream:
            stream.extractall(out / dest, filter="data")
    arms = {k: getattr(args, k).resolve() for k in ("parent", "linear", "cubic")}
    if args.shift:
        arms["shift"] = args.shift.resolve()
    orders = itertools.cycle(itertools.permutations(arms))
    rows = []
    for case in sorted((out / "starts").iterdir()):
        fixture = out / "fixtures" / case.name
        info = json.loads((fixture / "metadata.json").read_text())
        serial, refs = info["serial"], np.asarray(info["mapped"])
        common = [
            "--receptor",
            fixture / "receptor.pdbqt",
            "--seed",
            "63419",
            "--cpu",
            "1",
        ]
        for axis, value in zip("xyz", info["center"]):
            common.extend([f"--center_{axis}", str(value), f"--size_{axis}", "24"])
        for label in ("projected", "63419", "63421"):
            start = case / (label + ".pdbqt")
            start_rmsd = rms(read_pose(start, serial), refs)
            start_e = common_score(
                arms["parent"], common, fixture / "ligand.pdbqt", start
            )
            for cap in (12, 36):
                for arm in next(orders):
                    stem = f"{case.name}-{label}-{cap}-{arm}"
                    pose = out / "raw" / (stem + ".pdbqt")
                    pose.parent.mkdir(exist_ok=True)
                    before = perf_counter()
                    log = call(
                        [
                            arms[arm],
                            *common,
                            "--ligand",
                            start,
                            "--local_only",
                            "--refinement_optimizer",
                            "armijo",
                            "--max_evals",
                            str(cap),
                            "--out",
                            pose,
                        ]
                    )
                    seconds = perf_counter() - before
                    (out / "raw" / (stem + ".log")).write_text(log)
                    evaluations, objective = re.search(
                        r"gradient_evals=(\d+) objective_kcal_mol=(\S+)", log
                    ).groups()
                    if not (1 <= int(evaluations) <= cap):
                        raise ValueError("Violated local gradient budget: " + stem)
                    rows.append(
                        {
                            "pdb": case.name,
                            "target": info["target"],
                            "start": label,
                            "cap": cap,
                            "arm": arm,
                            "start_rmsd": start_rmsd,
                            "end_rmsd": rms(read_pose(pose, serial), refs),
                            "initial_common_energy": start_e,
                            "final_common_energy": common_score(
                                arms["parent"], common, fixture / "ligand.pdbqt", pose
                            ),
                            "gradient_evals": int(evaluations),
                            "objective": float(objective),
                            "seconds": round(seconds, 6),
                        }
                    )
    if len(rows) != 48 * len(arms):
        raise ValueError("Incomplete eight-target panel")
    with (out / "local.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
