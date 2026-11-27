"""Re-run all five ligands against both receptors with the final production binary."""

import argparse
import json
import shutil
from pathlib import Path

from assay import MEASURED, OLD, read, run_panel


def poses(path):
    models = []
    for block in path.read_text().split("ENDMDL")[:-1]:
        energy = next(
            line.split()[3]
            for line in block.splitlines()
            if line.startswith("REMARK VINA RESULT:")
        )
        atoms = tuple(
            line for line in block.splitlines() if line.startswith(("ATOM  ", "HETATM"))
        )
        models.append((energy, atoms))
    return models


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("work", "original", "new", "release"):
        p.add_argument("--" + key, type=Path, required=True)
    p.add_argument("--cpu", type=int, default=4)
    args = p.parse_args()
    args.reverse = False
    args.rescore = False
    # A final-release search is counted directly by replay of the same binary.
    # No upstream process is used in this single-arm verification.
    args.counter = args.release
    args.work.mkdir(exist_ok=True, parents=True)
    for ligand in MEASURED:
        code = ligand["pdb"]
        shutil.copytree(
            (args.original if code in OLD else args.new) / code,
            args.work / code,
            dirs_exist_ok=True,
        )
    run_panel.run(args, systems=[l["pdb"] for l in MEASURED], arms=("release",))
    prior = {
        r["key"]: r
        for path in (args.original, args.new)
        for r in read(path / "cells.csv")
        if r["arm"] == "candidate" and r["target"] == "FABP4"
    }
    checks = []
    for r in read(args.work / "cells.csv"):
        k = r["key"]
        pk = k.replace("_release", "_candidate")
        old = (
            (args.original if r["query"] in OLD else args.new)
            / "runs"
            / (pk + ".pdbqt")
        )
        now = args.work / "runs" / (k + ".pdbqt")
        a, b = poses(old), poses(now)
        assert len(a) == len(set(a)) and len(b) == len(set(b))
        same_work = all(
            prior[pk][i] == r[i]
            for i in (
                "hunt_calls",
                "hunt_evals",
                "second_calls",
                "second_evals",
                "explicit_calls",
                "explicit_evals",
                "explicit_inputs",
            )
        )
        checks.append(
            dict(
                key=k,
                identical_top1=a[0] == b[0],
                identical_work=same_work,
                identical_output=old.read_bytes() == now.read_bytes(),
                previous_modes=len(a),
                final_modes=len(b),
                removed_ranks=" ".join(
                    str(i) for i, m in enumerate(a, 1) if m not in b
                ),
                added_final_ranks=" ".join(
                    str(i) for i, m in enumerate(b, 1) if m not in a
                ),
            )
        )
    run_panel.paired.write(args.work / "comparisons.csv", checks)
    report = dict(
        release_sha256=run_panel.sha(args.release),
        cells=len(checks),
        same_top1=sum(r["identical_top1"] for r in checks),
        same_work=sum(r["identical_work"] for r in checks),
        same_bytes=sum(r["identical_output"] for r in checks),
        removed_modes=sum(len(r["removed_ranks"].split()) for r in checks),
        added_modes=sum(len(r["added_final_ranks"].split()) for r in checks),
    )
    (args.work / "release_verification.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
