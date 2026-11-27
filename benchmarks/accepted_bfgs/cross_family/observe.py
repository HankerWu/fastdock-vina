"""Verify native accepted energy and FG count on byte-identical instrumentation replays."""

import argparse
import csv
import json
import re
import tarfile
from pathlib import Path

from replay import digest, run, write_rows


def observe(work, binaries):
    raw = work / "trace"
    raw.mkdir()
    rows = []
    with (work / "rows.csv").open() as stream:
        for row in csv.DictReader(stream):
            arm, pdb, start = (row[k] for k in ("arm", "pdb", "start"))
            prefix = f"trace/{pdb}-{start}-{arm}"
            end = Path(prefix + ".pdbqt")
            log = work / (prefix + ".log")
            fixture = Path("fixtures") / pdb
            center = json.loads((work / fixture / "metadata.json").read_text())[
                "center"
            ]
            args = [
                "--receptor",
                fixture / "receptor.pdbqt",
                "--ligand_bonds",
                f"raw/{pdb}.bonds",
                "--cpu",
                "1",
                "--seed",
                "63419",
                "--ligand",
                f"starts/{pdb}/{start}.pdbqt",
                "--local_only",
                "--out",
                end,
            ]
            for axis, value in zip("xyz", center, strict=True):
                args.extend([f"--center_{axis}", str(value), f"--size_{axis}", "24"])
            run(binaries[arm], args, log, work)
            if digest(work / end) != row["end_sha256"]:
                raise ValueError(f"Instrumented native output differs in {prefix}")
            count, objective = re.findall(
                r"Local objective observation: gradient_evals=(\d+) accepted_objective=(\S+)",
                log.read_text(),
            )[0]
            rows.append(
                {
                    "pdb": pdb,
                    "start": start,
                    "arm": arm,
                    "native_fg": int(count),
                    "returned_native_objective": float(objective),
                    "verified_end_sha256": row["end_sha256"],
                }
            )
    write_rows(work / "native.csv", rows)
    (work / "observed_binaries.json").write_text(
        json.dumps({a: digest(b) for a, b in binaries.items()}, indent=2) + "\n"
    )
    with tarfile.open(work / "trace.tar.xz", "w:xz") as tar:
        tar.add(raw, arcname="trace")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "observed_parent", "observed_fixed"):
        ap.add_argument(name, type=lambda value: Path(value).resolve())
    args = ap.parse_args()
    observe(args.work, {"parent": args.observed_parent, "fixed": args.observed_fixed})
