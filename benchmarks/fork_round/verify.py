"""Rebuild this round's evidence from a clean checkout. Historical audits have a separate entry point."""
import filecmp
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
ROOT = BENCH.parent
os.environ["OPENBLAS_NUM_THREADS"] = "1"


def run(script, *args):
    subprocess.run([sys.executable, script, *args], cwd=ROOT, check=True)


def equal(actual, expected):
    if not filecmp.cmp(actual, expected, shallow=False):
        raise ValueError(f"Reconstruction differs: {expected}")


def verify(work):
    props = BENCH / "proposal_fork"
    merge = BENCH / "merge_replay"
    structure = BENCH / "fork_structure"
    assays = HERE / "affinity"
    for name, evidence in (("replay.py","raw"), ("rescue.py","rescue"), ("free.py","free")):
        run(props / name, props / evidence)
    run(props / "check.py")
    run(props / "gate.py")
    run(merge / "replay.py", "analyze", merge / "raw", work / "merge")
    reference_tables = {p.name for p in (merge/"results").glob("*.csv")}
    assert reference_tables == {p.name for p in (work/"merge").glob("*.csv")}
    for file in reference_tables:
        equal(work/"merge"/file, merge/"results"/file)
    run(merge/"verify.py", work/"merge.json")
    equal(work/"merge.json", merge/"verification.json")
    run(merge/"refinement.py", merge/"refinement")
    run(structure/"verify.py")
    run(BENCH/"search_trial_affinity/study.py", "replay", work/"assays", "--artifacts", assays)
    new = work/"assays/recomputed"
    run(BENCH/"stage_local_affinity/paired.py", new)
    run(BENCH/"stage_local_affinity/costs.py", work/"assays", new/"physical_costs.json")
    for file in ("pair_changes.csv","pair_summary.csv","bace1_receptor_pairs.csv","physical_costs.json"):
        equal(new/file, assays/file)
    with tarfile.open(HERE/"production/raw.tar.xz") as archive:
        archive.extractall(work/"production", filter="data")
    run(HERE/"production.py", work/"production", work/"assays", HERE/"production")
    sources = json.loads((structure/"builds.json").read_text())["versions"]
    for arm, study in json.loads((new/"physical_costs.json").read_text()).items():
        assert sources[arm]["binary_sha256"] == study["binary_sha256"]
    for arm in ("parent","stable","candidate"):
        rev = sources[arm]["revision"]
        tree = subprocess.check_output(["git","rev-parse",rev+":src"],cwd=ROOT,text=True).strip()
        assert tree == sources[arm]["src_tree"]
    run(BENCH/"cumulative/source/verify.py")
    run(HERE/"decision.py")
    run(BENCH/"cumulative/ledger.py")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as directory:
        verify(Path(directory))
    print("Verified locked and free paths, native merge, fixed and external poses, resource receipts, assay pairs and promotion.")
