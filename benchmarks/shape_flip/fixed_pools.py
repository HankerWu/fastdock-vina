"""Check all three Pro builds on the complete unmodified development-parent pools."""

import argparse
import csv
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

EXTERNAL = Path(__file__).resolve().parents[1] / "macrocycle_proposal" / "external"
sys.path.insert(0, str(EXTERNAL))
from external_run import write


def main(a):
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        pose = Path(tmp) / "pose.pdbqt"
        for r in csv.DictReader((a.base / "runs.csv").open()):
            if r["arm"] != "parent":
                continue
            directory = a.inputs / r["pdb"]
            center = json.loads((directory / "preparation.json").read_text())["center"]
            box = [
                v
                for i, c in zip("xyz", center)
                for v in ("--center_" + i, str(c), "--size_" + i, "26")
            ]
            blocks = (
                (a.base / "raw" / (r["key"] + ".pdbqt"))
                .read_text()
                .split("ENDMDL")[:-1]
            )
            for rank, block in enumerate(blocks, 1):
                pose.write_text(block.strip() + "\nENDMDL\n")
                terms = []
                for binary in (a.parent, a.always, a.adaptive):
                    result = subprocess.run(
                        [
                            binary,
                            "--score_only",
                            "--receptor",
                            directory / "receptor.pdbqt",
                            "--ligand",
                            directory / "ligand.pdbqt",
                            "--pose_coordinates",
                            pose,
                            "--unbound_energy",
                            "0",
                            "--seed",
                            r["seed"],
                            "--cpu",
                            "1",
                            *box,
                        ],
                        capture_output=True,
                        check=False,
                        text=True,
                    )
                    if result.returncode:
                        if "outside the grid box" not in result.stderr:
                            raise ValueError(result.stderr)
                        terms.append(["outside"])
                    else:
                        values = re.findall(
                            r":\s*([-+.\deE]+)\s+\(kcal/mol\)", result.stdout
                        )
                        if len(values) != 10:
                            raise ValueError(result.stdout)
                        terms.append(values)
                if not terms[0] == terms[1] == terms[2]:
                    raise ValueError(r["key"])
                rows.append(
                    {"key": r["key"], "rank": rank, "terms": " ".join(terms[0])}
                )
    write(a.output, rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for s in ("inputs", "base", "output", "parent", "always", "adaptive"):
        p.add_argument("--" + s, type=Path, required=True)
    main(p.parse_args())
