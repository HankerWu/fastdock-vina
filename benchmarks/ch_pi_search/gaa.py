"""Prepare the two GAA structures and execute the prespecified cross-docks."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
from study import run_cases
from experiment import prepare, cross

CASES = {"5NN5": ("A", "NOJ"), "5NN6": ("A", "MIG")}


def dock(work, binary, out):
    meta = json.loads((work / "metadata.json").read_text())
    jobs = []
    alignments = []
    for src, dst in product(CASES, repeat=2):
        rot, shift, residual, count = cross.alignment(meta[src], meta[dst])
        refs = np.asarray(meta[src]["mapped"]) @ rot + shift
        alignments.append(
            dict(query=src, receptor=dst, pocket_ca=count, residual=residual)
        )
        for seed, cap in product(cross.SEEDS, (12000, 24000)):
            args = [
                "--ligand",
                f"{src}/protonated.pdbqt",
                "--receptor",
                f"{dst}/receptor.pdbqt",
                "--cpu",
                "1",
                "--seed",
                str(seed),
                "--exhaustiveness",
                "4",
                "--max_evals",
                str(cap),
                "--num_modes",
                "9",
                "--energy_range",
                "100",
                "--pose_selection",
                "late",
                "--refine_candidates",
                "40",
                "--search_counters",
                "--out",
                "pose.pdbqt",
            ]
            for axis, center in zip("xyz", meta[dst]["center"]):
                args.extend(["--center_" + axis, str(center), "--size_" + axis, "24"])
            jobs.append(
                (
                    ("GAA", src, dst, seed, cap),
                    work,
                    args,
                    None,
                    meta[src]["protonated_serial"],
                    refs,
                )
            )
    run_cases(jobs, binary, out)
    (out / "alignment.json").write_text(json.dumps(alignments, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=lambda s: Path(s).resolve())
    sub = p.add_subparsers(dest="action", required=True)
    sub.add_parser("prepare")
    dock_parser = sub.add_parser("dock")
    for arg in ("binary", "out"):
        dock_parser.add_argument(arg, type=lambda s: Path(s).resolve())
    a = p.parse_args()
    if a.action == "prepare":
        prepare(a.work, CASES)
    else:
        dock(a.work, a.binary, a.out)
