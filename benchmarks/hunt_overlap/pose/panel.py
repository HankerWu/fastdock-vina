"""Prospective THRB/BRD9 pose experiment using the repository's paired evaluator."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
import sys
from pathlib import Path

CASES = {
    "1NAX": ("THRB", "IH5", "A"),
    "3IMY": ("THRB", "B72", "A"),
    "5EU1": ("BRD9", "5SW", "A"),
    "5F1H": ("BRD9", "5U6", "A"),
}
HOLDOUT = {"1UPV": ("NR1H2", "444", "A"), "1P8D": ("NR1H2", "CO1", "A")}
CONFIRM = {**HOLDOUT, "3VJK": ("DPP4", "M51", "A"), "4A5S": ("DPP4", "N7F", "A")}
STUDIES = ("tangent", "placement-dev", "placement-holdout", "placement-confirm")


def design(study):
    if study == "placement-confirm":
        return CONFIRM, (704521,), ("upstream", "parent", "fit", "invalid")
    if study == "tangent":
        return CASES, (481211, 481217), ("upstream", "parent", "all", "first")
    cases, seeds = (
        (HOLDOUT, (625231, 625237))
        if study == "placement-holdout"
        else (CASES, (481211, 481217))
    )
    return cases, seeds, ("upstream", "parent", "fit", "invalid")


def configure_arm(arm):
    os.environ.pop("FASTDOCK_HUNT_TANGENT", None)
    os.environ.pop("FASTDOCK_FIRST_PLACEMENT", None)
    if arm in ("all", "first"):
        os.environ["FASTDOCK_HUNT_TANGENT"] = arm
    if arm in ("fit", "invalid"):
        os.environ["FASTDOCK_FIRST_PLACEMENT"] = arm


CAPS = (3000, 12000)


def graph_and_receptor_tools(repo, cases):
    sys.path.insert(0, str(repo / "benchmarks/search_backtracking/pose"))
    import backtracking_panel as panel
    import experiment

    panel.CASES = experiment.CASES = cases
    return panel, experiment


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "run", "verify", "parity", "crystal"))
    p.add_argument("work", type=Path)
    p.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3])
    p.add_argument("--study", choices=STUDIES, default="tangent")
    p.add_argument("--core", type=Path)
    p.add_argument("--upstream", type=Path)
    p.add_argument("--counts", type=Path)
    p.add_argument("--release", type=Path)
    p.add_argument("--pristine", type=Path)
    p.add_argument("--cpu", type=int)
    a = p.parse_args()
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "OPENMM_CPU_THREADS"):
        os.environ[key] = "1"
    if a.cpu is not None:
        os.sched_setaffinity(0, {a.cpu})
    cases, seeds, arms = design(a.study)
    panel, experiment = graph_and_receptor_tools(a.repo.resolve(), cases)
    paired = panel.paired
    work = a.work.resolve()
    if a.action == "prepare":
        for pdb in cases:
            panel.prepare_input(work, pdb)
        panel.write_sources(work)
        return
    if a.action == "verify":
        import audit

        result = audit.main(work, cases)
        (work / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
        return
    meta = {q: json.loads((work / q / "preparation.json").read_text()) for q in cases}
    pairs = experiment.geometry(work, meta)
    if a.action == "crystal":
        if not a.pristine:
            p.error("crystal requires --pristine")
        import numpy as np
        from crystal_start import run as relax

        (work / "metadata.json").write_text(json.dumps(meta) + "\n")
        # --local_only requires the input ligand to be inside the identical search box.
        # Preserve every blind cell, including receivers with an out-of-box reference.
        starts = []
        eligible = set()
        for q, r in pairs:
            rot, shift, _, _ = paired.base.alignment(meta[q], meta[r])
            delta = np.asarray(meta[q]["mapped"][0]) @ rot + shift - meta[r]["center"]
            widths = abs(delta).max(axis=1)
            outside = int((widths > 12).sum())
            starts.append(
                {
                    "query": q,
                    "receiver": r,
                    "center_distance_A": float(np.linalg.norm(delta.mean(axis=0))),
                    "max_extent_A": float(widths.max()),
                    "outside_heavy_atoms": outside,
                }
            )
            if outside == 0:
                eligible.add((q, r))
        paired.write(work / "crystal_box.csv", starts)
        groups = {False: set(), True: set()}
        for q, r in eligible:
            groups[meta[q]["target"] != "NR1H2"].add((q, r))
        outputs = []
        for compare_mappings, group in groups.items():
            if group:
                relax(
                    work, a.pristine.resolve(), group, compare_mappings=compare_mappings
                )
                with (work / "crystal_starts.csv").open() as handle:
                    outputs.extend(csv.DictReader(handle))
        paired.write(work / "crystal_starts.csv", outputs)
        return
    cells = [(q, r, s, c) for q, r in pairs for s, c in itertools.product(seeds, CAPS)]
    random.Random(seeds[0]).shuffle(cells)
    (work / "runs").mkdir(exist_ok=True)
    if a.action == "parity":
        if not a.pristine or not a.release:
            p.error("parity requires --pristine and --release")
        controls = []
        for q, r, s, c in cells:
            for original, arm, binary in (
                ("upstream", "release", a.release),
                ("parent", "pristine", a.pristine),
            ):
                configure_arm(arm)
                row, _ = paired.execute(work, meta, q, r, s, c, arm, binary, False)
                reference = work / "runs" / f"{q}_{r}_{s}_{c}_{original}.pdbqt"
                same = digest(reference) == digest(
                    work / "runs" / (row["key"] + ".pdbqt")
                )
                if not same:
                    raise ValueError("Byte parity failed for " + row["key"])
                controls.append(
                    {
                        "key": row["key"],
                        "comparator": original,
                        "identical": same,
                        "pdbqt_sha256": digest(reference),
                        "binary_sha256": digest(binary),
                        "seconds": row["seconds"],
                    }
                )
        paired.write(work / "parity.csv", controls)
        return
    if not all((a.core, a.upstream, a.counts)):
        p.error("run requires --core, --upstream and --counts")
    rows, modes = [], []
    for i, (q, r, s, c) in enumerate(cells):
        order = i % len(arms)
        for arm in arms[order:] + arms[:order]:
            configure_arm(arm)
            binary = a.upstream if arm == "upstream" else a.core
            row, poses = paired.execute(
                work, meta, q, r, s, c, arm, binary, False, a.counts
            )
            for field in (
                "hunt_evals",
                "second_evals",
                "explicit_evals",
                "hunt_calls",
                "second_calls",
                "explicit_calls",
            ):
                if row[field] == "":
                    raise ValueError(f"Missing {field} in {row['key']}")
            row["pdbqt_sha256"] = digest(work / "runs" / (row["key"] + ".pdbqt"))
            rows.append(row)
            modes.extend(poses)
    paired.write(work / "cells.csv", rows)
    paired.write(work / "modes.csv", modes)
    ids = {k: digest(getattr(a, k)) for k in ("core", "counts", "upstream")}
    (work / "binaries.json").write_text(json.dumps(ids, indent=2) + "\n")


if __name__ == "__main__":
    main()
