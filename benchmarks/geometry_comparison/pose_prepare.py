"""Prepare the frozen human structures with the existing geometry-panel protocol."""

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "geometry_metric"))
import panel


def prepare(job):
    work, case = job
    case = dict(case)
    pdb = case["pdb"]
    out = work / pdb
    try:
        meta = panel.prepare(work, pdb, case["ligand"], case["target"])
        for filename in ("ligand.pdbqt", "receptor.pdbqt", "metadata.json"):
            case[filename + "_sha256"] = hashlib.sha256(
                (out / filename).read_bytes()
            ).hexdigest()
        case["status"] = meta["status"]
    except (ValueError, RuntimeError, KeyError, subprocess.CalledProcessError) as exc:
        case["status"] = "excluded"
        case["reason"] = str(exc).replace(str(work), "WORK")
    for filename in (pdb + ".pdb", case["ligand"] + ".cif"):
        if (out / filename).is_file():
            case[filename + "_sha256"] = hashlib.sha256(
                (out / filename).read_bytes()
            ).hexdigest()
    if (out / "source.json").is_file():
        case["source"] = json.loads((out / "source.json").read_text())
    print(pdb, case["status"], case.get("reason", ""), flush=True)
    return case


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    args = p.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    with (HERE / "pose_cases.csv").open() as f:
        cases = list(csv.DictReader(f))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = list(pool.map(prepare, ((args.work.resolve(), case) for case in cases)))
    (args.work / "manifest.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
