"""Score the identical Gal3 top poses and the complete structural development output pools."""

import argparse
import csv
import json
import re
import subprocess
import tempfile
from pathlib import Path


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main(args):
    jobs = []
    for row in read(args.structure / "cells.csv"):
        if row["arm"] != "parent":
            continue
        rec = row["receiver"]
        site = json.loads((args.structure / rec / "preparation.json").read_text())
        flags = ["--receptor", str(args.structure / rec / "receptor.pdbqt")]
        for axis, center in zip("xyz", site["center"]):
            flags.extend(["--center_" + axis, str(center), "--size_" + axis, "24"])
        jobs.append(
            (
                row["key"],
                args.structure / "runs" / (row["key"] + ".pdbqt"),
                flags,
                False,
            )
        )
    for row in read(args.affinity / "cells.csv"):
        if row["arm"] != "parent":
            continue
        key = row["receptor"] + "_" + row["compound"] + "_parent"
        rec = row["receptor"]
        flags = [
            "--receptor",
            str(args.inputs / (rec + ".pdbqt")),
            "--config",
            str(args.inputs / (rec + ".box")),
        ]
        jobs.append((key, args.affinity / "runs" / (key + ".pdbqt"), flags, True))
    with tempfile.TemporaryDirectory() as tmp, args.csv.open("w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(
            [
                "case",
                "mode",
                "status",
                "fixed_reference_score",
                "inter",
                "intra",
                "torsions",
                "unbound",
            ]
        )
        path = Path(tmp) / "fixed.pdbqt"
        for key, pool, flags, first in jobs:
            blocks = pool.read_text().split("ENDMDL")[:-1]
            for i, block in enumerate(blocks[:1] if first else blocks, 1):
                path.write_text(
                    "\n".join(
                        s for s in block.splitlines() if not s.startswith("MODEL")
                    )
                    + "\n"
                )
                terms = []
                for executable in (args.parent, args.always, args.adaptive):
                    run = subprocess.run(
                        [
                            str(executable),
                            *flags,
                            "--cpu",
                            "1",
                            "--seed",
                            "214039",
                            "--ligand",
                            str(path),
                            "--score_only",
                            "--unbound_energy",
                            "0",
                        ],
                        text=True,
                        capture_output=True,
                    )
                    if run.returncode:
                        if "outside the grid box" not in run.stderr + run.stdout:
                            raise RuntimeError((key, i, run.stderr))
                        terms.append(("outside",))
                        continue
                    output = run.stdout.split("Estimated Free Energy of Binding")[-1]
                    energies = re.findall(r":\s+(-?\d+\.\d+)\s+\(kcal/mol\)", output)
                    if len(energies) != 10:
                        raise ValueError(output)
                    terms.append(("ok", *energies))
                assert terms[0] == terms[1] == terms[2], (key, i, terms)
                result = terms[0]
                selected = (
                    result
                    if result[0] == "outside"
                    else [result[k] for k in (0, 1, 2, 5, 9, 10)]
                )
                writer.writerow([key, i, *selected])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "structure",
        "affinity",
        "inputs",
        "parent",
        "always",
        "adaptive",
        "csv",
    ):
        parser.add_argument(name, type=Path)
    main(parser.parse_args())
