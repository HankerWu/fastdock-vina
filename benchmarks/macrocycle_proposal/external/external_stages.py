"""Replay the sampled, proposed, minimized and retained poses; verify observer parity."""

import argparse
import csv
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path

from external_run import sha256, write
from external_score import main as score


def main(work, base, out, binaries):
    raw = out / "raw"
    raw.mkdir(exist_ok=True, parents=True)
    records = []
    for r in csv.DictReader((base / "runs.csv").open()):
        if r["arm"] not in binaries or r["cap"] != "3000":
            continue
        d = work / r["pdb"]
        center = json.loads((d / "preparation.json").read_text())["center"]
        box = [
            x
            for a, c in zip("xyz", center)
            for x in ("--center_" + a, str(c), "--size_" + a, "26")
        ]
        prefix = raw / r["key"]
        if list(raw.glob(r["key"] + ".*.pdbqt")):
            raise ValueError("Existing stage outputs for " + r["key"])
        with tempfile.TemporaryDirectory() as tmp:
            replay = Path(tmp) / "out.pdbqt"
            args = [
                binaries[r["arm"]],
                "--receptor",
                d / "receptor.pdbqt",
                "--ligand",
                d / "ligand.pdbqt",
                *box,
                "--cpu",
                "1",
                "--seed",
                r["seed"],
                "--exhaustiveness",
                "1",
                "--max_evals",
                r["cap"],
                "--verbosity",
                "0",
                "--num_modes",
                "5",
                "--energy_range",
                "100",
                "--out",
                replay,
                "--search_counters",
            ]
            log = subprocess.check_output(
                args,
                text=True,
                stderr=subprocess.STDOUT,
                env={**os.environ, "FASTDOCK_RING_TRACE": str(prefix.resolve())},
            )
            (raw / (r["key"] + ".log")).write_text(log)
            counts = dict(
                re.findall(
                    r"(hunt_calls|hunt_evals|second_calls|second_evals|explicit_calls|explicit_evals)=(\d+)",
                    log,
                )
            )
            if (
                len(counts) != 6
                or any(r[k] != v for k, v in counts.items())
                or sha256(replay) != r["sha256"]
            ):
                raise ValueError("Observer alters output " + r["key"])
        for stage in ("raw", "proposal", "hunt", "second", "saved"):
            key = r["key"] + "." + stage
            n = (raw / (key + ".pdbqt")).read_text().count("ENDMDL")
            expected = int(r.get(stage + "_calls", 1))
            if n == 0 or (stage != "saved" and n != expected):
                raise ValueError("Invalid stage count " + key)
            records.append(
                {
                    "key": key,
                    "pdb": r["pdb"],
                    "arm": r["arm"],
                    "cap": r["cap"],
                    "seed": r["seed"],
                    "stage": stage,
                    "modes": n,
                }
            )
    write(out / "runs.csv", records)
    score(work, out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "base", "out"):
        p.add_argument(k, type=Path)
    for k in ("parent", "reset", "ccd", "lock"):
        p.add_argument("--" + k, type=Path)
    a = p.parse_args()
    binaries = {
        k: getattr(a, k)
        for k in ("parent", "reset", "ccd", "lock")
        if getattr(a, k) is not None
    }
    if not binaries:
        p.error("Select at least one observer binary")
    main(a.work, a.base, a.out, binaries)
