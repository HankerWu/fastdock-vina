"""Run native paths, replay both proposal libraries through both historical binaries."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PANEL = HERE.parent / "geometry_comparison/upstream/four_versions"
sys.path.insert(0, str(PANEL))
import run_panel

ARMS = ("e7031e9", "2326430")
CASES = {("2NNQ", "3FR2", 85237), ("1FKG", "1FKH", 85237), ("5TBM", "4GHI", 85231)}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(bins, work):
    work.mkdir(parents=True, exist_ok=False)
    costs = []
    archive = list(csv.DictReader((PANEL / "bfgs_results/all_cells.csv").open()))
    for job in run_panel.panel(work / "fixtures"):
        cohort, inputs, q, r, reference, seed, cap = job
        case = (q["pdb"], r["pdb"], seed)
        if cohort != "reciprocal" or case not in CASES:
            continue
        path = work / ("-".join(map(str, case)))
        path.mkdir()
        for arm, binary in bins.items():
            os.environ.pop("FD_BFGS_SHADOW", None)
            os.environ["FD_BFGS_TRACE"] = str(path / (arm + ".native"))
            stats, ranks = run_panel.dock((arm, binary, path, job))
            assert stats["status"] == "ok"
            historical = next(row for row in archive if row["key"] == stats["key"])
            assert stats["pool_sha256"] == historical["pool_sha256"]
            native = Path(os.environ["FD_BFGS_TRACE"])
            proposals = path / (arm + ".inputs")
            with proposals.open("w") as dest:
                for line in native.read_text().splitlines():
                    if line.startswith("P\t"):
                        print(line.split("\t")[3], file=dest)
            costs.append(dict(case=path.name, driver=arm, evaluator=arm, kind="native",
                              trace_sha256=digest(native), binary_sha256=digest(binary), **stats))
        for arm in ARMS:
            for evaluator, binary in bins.items():
                # Panel executable, receptor, grid and hunt capacity are identical to
                # the free runs. Observer exits after solving the supplied inputs.
                prefix = f"{arm}-{evaluator}"
                replay = path / (prefix + ".shadow")
                os.environ["FD_BFGS_SHADOW"] = str(path / (arm + ".inputs"))
                os.environ["FD_BFGS_TRACE"] = str(replay)
                cmd = [str(binary), "--ligand", str(inputs / q["pdb"] / "ligand.pdbqt"),
                       "--receptor", str(inputs / r["pdb"] / "receptor.pdbqt"), "--seed", str(seed),
                       "--cpu", "1", "--max_evals", str(cap), "--exhaustiveness", "3",
                       "--num_modes", "9", "--energy_range", "100", "--spacing", "0.375",
                       "--out", str(path / (prefix + ".unused.pdbqt"))]
                for axis, value in zip("xyz", r["center"]):
                    cmd += ["--size_" + axis, "24", "--center_" + axis, str(value)]
                stats = run_panel.execute(cmd, path / (prefix + ".log"), 180)
                assert stats["status"] == "ok"
                costs.append(dict(case=path.name, driver=arm, evaluator=evaluator, kind="shadow",
                                  trace_sha256=digest(replay), binary_sha256=digest(binary), **stats))
    (work / "costs.json").write_text(json.dumps(costs, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("old", type=Path)
    parser.add_argument("repaired", type=Path)
    args = parser.parse_args()
    run(dict(zip(ARMS, (args.old.resolve(), args.repaired.resolve()))), args.output.resolve())
