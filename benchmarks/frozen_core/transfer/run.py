"""Prepare two human target pairs and generate one native pool per query, receiver and seed."""

import argparse
import hashlib
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "geometry_metric"))
import panel

CASES = (
    ("PARP1", "7KK4", "09L"),
    ("PARP1", "4UND", "2YQ"),
    ("CDK2", "2A4L", "RRC"),
    ("CDK2", "1AQ1", "STU"),
)
SEEDS = (98207, 98213)


def prepare(work):
    report = []
    for target, pdb, ligand in CASES:
        try:
            example = panel.prepare(work, pdb, ligand, target)
            report.append(
                dict(target=target, pdb=pdb, ligand=ligand, status=example["status"])
            )
        except Exception as error:
            report.append(
                dict(
                    target=target,
                    pdb=pdb,
                    ligand=ligand,
                    status="excluded",
                    reason=str(error),
                )
            )
        print(report[-1], flush=True)
    (work / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")


def dock(work, binary):
    out = work / "pools"
    out.mkdir(exist_ok=True)
    manifest = json.loads((work / "manifest.json").read_text())

    def one(job):
        receiver, query, seed = job
        rm = json.loads((work / receiver / "metadata.json").read_text())
        stem = out / f"{query}-{receiver}-{seed}"
        command = [
            str(binary.resolve()),
            "--ligand",
            str(work / query / "ligand.pdbqt"),
            "--receptor",
            str(work / receiver / "receptor.pdbqt"),
            "--seed",
            str(seed),
            "--cpu",
            "1",
            "--exhaustiveness",
            "4",
            "--max_evals",
            "60000",
            "--num_modes",
            "9",
            "--energy_range",
            "100",
            "--out",
            str(stem.with_suffix(".pdbqt")),
        ]
        for axis, center in zip("xyz", rm["center"]):
            command.extend(("--size_" + axis, "24", "--center_" + axis, str(center)))
        with stem.with_suffix(".log").open("w") as stream:
            subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)
        assert stem.with_suffix(".pdbqt").read_bytes().endswith(b"ENDMDL\n")
        return stem.name

    jobs = [
        (a["pdb"], b["pdb"], s)
        for a in manifest
        for b in manifest
        for s in SEEDS
        if a["status"] == b["status"] == "ready"
        and a["pdb"] != b["pdb"]
        and a["target"] == b["target"]
    ]
    with ThreadPoolExecutor(4) as jobs_pool:
        for name in jobs_pool.map(one, jobs):
            print(name, flush=True)
    (work / "binary.sha256").write_text(
        hashlib.sha256(binary.read_bytes()).hexdigest() + "  vina\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "dock"))
    parser.add_argument("work", type=Path)
    parser.add_argument("--binary", type=Path)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    if args.action == "prepare":
        prepare(args.work)
    else:
        dock(args.work, args.binary)
