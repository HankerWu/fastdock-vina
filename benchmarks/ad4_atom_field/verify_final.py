"""Rebuild confidence in the release: replay its 96 actual assay outputs and counters."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile
from tempfile import TemporaryDirectory


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("work", "binary", "out"):
        p.add_argument(arg, type=Path)
    p.add_argument("--jobs", type=int, default=2)
    a = p.parse_args()
    here = Path(__file__).resolve().parent
    tasks = []
    for cohort in ("affinity_primary", "affinity_confirmation"):
        src = here / cohort
        with tarfile.open(src / "pools.tar.xz") as archive:
            for row in csv.DictReader((src / "affinity_search.csv").open()):
                if row["arm"] != "release":
                    continue
                name = "-".join(row[k] for k in ("ligand", "receiver", "seed", "arm"))
                data = archive.extractfile("pools/" + name + ".pdbqt").read()
                tasks.append((row, hashlib.sha256(data).hexdigest()))

    def check(task):
        row, expected = task
        with TemporaryDirectory() as tmp:
            out = Path(tmp) / "out.pdbqt"
            r = subprocess.run(
                [
                    a.binary,
                    "--ligand",
                    a.work / "affinity" / (row["ligand"] + ".pdbqt"),
                    "--maps",
                    a.work / "maps" / row["receiver"] / "receptor",
                    "--scoring",
                    "ad4",
                    "--cpu",
                    "1",
                    "--exhaustiveness",
                    "3",
                    "--max_evals",
                    "24000",
                    "--seed",
                    row["seed"],
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--search_counters",
                    "--out",
                    out,
                ],
                text=True,
                capture_output=True,
                check=True,
            )
            counted = sum(
                int(t)
                for pair in re.findall(
                    r"hunt_evals=(\d+).*?second_evals=(\d+)", r.stdout
                )
                for t in pair
            )
            if hashlib.sha256(
                out.read_bytes()
            ).hexdigest() != expected or counted != int(row["evals"]):
                raise ValueError(row)
            return int(row["modes"]), counted

    with ThreadPoolExecutor(a.jobs) as pool:
        result = list(pool.map(check, tasks))
    a.out.write_text(
        json.dumps(
            dict(
                bytes_and_counts_equal=len(result),
                models=sum(m for m, e in result),
                evals=sum(e for m, e in result),
                sha256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
