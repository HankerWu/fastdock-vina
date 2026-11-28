"""Rebuild this round's evidence in a temporary directory from a clean checkout."""

import csv
import filecmp
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
ROOT = BENCH.parent


def run(script, *arguments):
    subprocess.run(
        [sys.executable, script, *arguments],
        cwd=ROOT,
        env={**os.environ, "OPENBLAS_NUM_THREADS": "1"},
        check=True,
    )


def compare(actual, expected):
    if not filecmp.cmp(actual, expected, shallow=False):
        raise ValueError(f"reconstruction differs: {expected}")


def unpack(source, dest):
    with tarfile.open(source) as archive:
        archive.extractall(dest, filter="data")


def restored_structure(work, structure):
    restored = structure / "restored"
    unpack(restored / "raw.tar.xz", work)

    def key(row):
        return tuple(row[k] for k in ("cohort", "query", "receiver", "seed"))

    initial = BENCH / "geometry_comparison/upstream/four_versions/results/cells.csv"
    with initial.open() as stream:
        originals = {
            key(r): r for r in csv.DictReader(stream) if r["arm"] == "candidate"
        }
    with (restored / "cells.csv").open() as stream:
        cells = list(csv.DictReader(stream))
    assert len(cells) == len(originals) == 65 and set(map(key, cells)) == set(originals)
    for row in cells:
        old = originals[key(row)]
        pose = (work / "raw" / (row["key"] + ".pdbqt")).read_bytes()
        assert (
            hashlib.sha256(pose).hexdigest() == row["pool_sha256"] == old["pool_sha256"]
        )
        assert all(
            row[k] == old[k]
            for k in (
                "modes",
                "status",
                "hunt_evals",
                "second_evals",
                "explicit_evals",
                "explicit_calls",
                "late_prescore",
            )
        )
    print(
        "Verified 65 restored production pools and all recorded stage counts",
        flush=True,
    )


def verify(work):
    hist = BENCH / "bfgs_source_replay"
    stage = BENCH / "stage_local_attribution"
    structure = BENCH / "stage_local_structure"
    assay = BENCH / "stage_local_affinity"
    study = BENCH / "search_trial_affinity/study.py"
    run(BENCH / "geometry_comparison/upstream/four_versions/offline.py")
    run(
        BENCH / "search_trial_panel/panel.py",
        "check",
        BENCH / "search_trial_panel/results",
    )
    run(hist / "archive.py", "restore", hist / "raw", work / "history")
    run(hist / "analyze.py", work / "history", work / "old-tables")
    run(hist / "candidate_replay.py", work / "history", work / "old-tables")
    for expected in (hist / "results").glob("*.csv"):
        compare(work / "old-tables" / expected.name, expected)
    run(hist / "verify.py", work / "history", hist / "results")

    run(stage / "pack.py", "restore", stage / "raw", work / "paths")
    run(stage / "analyze_stages.py", work / "paths", work / "tables")
    for expected in (stage / "results").glob("*.csv"):
        compare(work / "tables" / expected.name, expected)
    run(stage / "verify_stage.py", work / "paths", work / "tables")
    run(HERE / "locked.py", work / "history", HERE / "locked")
    run(HERE / "coverage.py", HERE / "coverage")

    run(structure / "cost.py", "check", structure / "fixed")
    restored_structure(work / "restored-structure", structure)
    unpack(structure / "external/inputs.tar.xz", work / "esr1")
    run(structure / "external.py", "check", work / "esr1", structure / "external")
    run(study, "replay", work / "assays", "--artifacts", assay / "results")
    recalculated = work / "assays/recomputed"
    run(assay / "paired.py", recalculated)
    run(assay / "costs.py", work / "assays", recalculated / "physical_costs.json")
    for name in (
        "pair_changes.csv",
        "pair_summary.csv",
        "bace1_receptor_pairs.csv",
        "physical_costs.json",
    ):
        compare(recalculated / name, assay / "results" / name)
    unpack(assay / "restored/restored.tar.xz", work / "production")
    run(
        assay / "restored.py",
        work / "production",
        work / "assays",
        work / "restored-check",
    )
    compare(work / "restored-check/restored.json", assay / "restored/restored.json")
    run(BENCH / "cumulative/source/verify.py")
    run(BENCH / "cumulative/ledger.py")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as directory:
        verify(Path(directory))
    print(
        "Verified historical and current four-source evidence, locked and free paths, transfer, assays and ledger."
    )
