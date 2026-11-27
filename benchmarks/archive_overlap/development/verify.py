"""Re-evaluate archived coordinates, work and full-precision candidate streams."""

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

from dock_panel import BENCH, counts, holdout
import replay

HERE = Path(__file__).resolve().parent


def extract(path, dest):
    with tarfile.open(path) as archive:
        archive.extractall(dest, filter="data")


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def verify(args):
    with tempfile.TemporaryDirectory() as folder:
        root = Path(folder)
        extract(BENCH / "rigid_flip/inputs.tar.xz", root)
        meta = {
            q: json.loads((root / q / "preparation.json").read_text())
            for q in holdout.CASES
        }
        total = nmode = rebuilt = 0
        for panel in (HERE, HERE / "factorial", HERE.parent / "refill"):
            work = root / panel.name
            work.mkdir()
            for q in meta:
                (work / q).symlink_to(root / q, target_is_directory=True)
            extract(panel / "raw.tar.xz", work)
            expected_modes = {
                (r["key"], int(r["mode"])): r for r in rows(panel / "modes.csv")
            }
            for row in rows(panel / "cells.csv"):
                key = row["key"]
                prefix = work / "raw" / key
                q, r = row["query"], row["receiver"]
                data = prefix.with_suffix(".pdbqt").read_bytes()
                if hashlib.sha256(data).hexdigest() != row["sha256"]:
                    raise ValueError(key + " output digest")
                logged = counts(prefix.with_suffix(".log").read_text())
                if any(int(row[k]) != v for k, v in logged.items()):
                    raise ValueError(key + " work")
                got = holdout.paired.inspect(work / "raw", meta[q], meta[r], dict(row))
                for pose in got:
                    wanted = expected_modes.pop((key, pose["mode"]))
                    if any(
                        abs(float(wanted[k]) - pose[k]) > 1e-9
                        for k in (
                            "rmsd_A",
                            "score",
                            "clash_pairs",
                            "minimum_heavy_distance_A",
                        )
                    ):
                        raise ValueError(key + " pose")
                total += 1
                nmode += len(got)
                if args.binary and row["arm"] == "fixed" and panel == HERE:
                    out = work / "replay.pdbqt"
                    cli = holdout.args_for(
                        root, meta, q, r, row["seed"], row["cap"], out
                    )
                    cli[cli.index("--exhaustiveness") + 1] = "4"
                    actual = subprocess.run(
                        [str(args.binary.resolve()), *cli, "--search_counters"],
                        capture_output=True,
                        text=True,
                        check=True,
                    )
                    if data != out.read_bytes() or counts(actual.stdout) != logged:
                        raise ValueError(key + " clean binary replay")
                    rebuilt += 1
            if expected_modes:
                raise ValueError("Extra archived modes")
            if panel == HERE:
                shutil.copyfile(panel / "cells.csv", work / "cells.csv")
                replay.main(argparse.Namespace(inputs=root, results=work))
                for name in ("streams.csv", "interventions.csv"):
                    if rows(work / name) != rows(panel / name):
                        raise ValueError("Frozen-candidate replay differs: " + name)
        print(
            f"{total} CLI outputs, {nmode} poses, {rebuilt} clean-build redocks verified"
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path)
    verify(p.parse_args())
