"""Score all parent poses of the topology-affected drug without moving coordinates."""

import argparse
import csv
import hashlib
import json
import re
import subprocess
import tempfile
import time
from pathlib import Path


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "parent", "fixed"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    w = a.work.resolve()
    mode_rows = {(r["key"], int(r["mode"])): r for r in rows(w / "modes.csv")}
    scores = []
    with tempfile.TemporaryDirectory() as temp:
        pose = Path(temp) / "pose.pdbqt"
        for run in rows(w / "cells.csv"):
            if run["arm"] != "parent":
                continue
            q, r, key = run["query"], run["receiver"], run["key"]
            source = json.loads((w / q / "preparation.json").read_text())
            if source["ligand"] != "3JW":
                continue
            receiver = json.loads((w / r / "preparation.json").read_text())
            block = (w / "runs" / (key + ".pdbqt")).read_text()
            for i, data in enumerate(
                re.findall(r"^MODEL .*?^ENDMDL", block, re.MULTILINE | re.DOTALL), 1
            ):
                pose.write_text(data + "\n")
                row = {
                    "key": key,
                    "query": q,
                    "receiver": r,
                    "cap": run["cap"],
                    "mode": i,
                    "rmsd_A": mode_rows[key, i]["rmsd_A"],
                    "pose_sha256": hashlib.sha256(pose.read_bytes()).hexdigest(),
                }
                for arm in ["parent", "fixed"]:
                    cmd = [
                        str(getattr(a, arm).resolve()),
                        "--score_only",
                        "--pose_coordinates",
                        str(pose),
                        "--unbound_energy",
                        "0",
                        "--cpu",
                        "1",
                        "--seed",
                        run["seed"],
                        "--ligand",
                        str(w / q / "ligand.pdbqt"),
                        "--receptor",
                        str(w / r / "receptor.pdbqt"),
                    ]
                    for axis, value in zip("xyz", receiver["center"]):
                        cmd += ["--center_" + axis, str(value), "--size_" + axis, "24"]
                    start = time.perf_counter()
                    log = subprocess.check_output(
                        cmd, text=True, stderr=subprocess.STDOUT
                    )
                    row[arm + "_seconds"] = time.perf_counter() - start
                    row[arm + "_kcal_mol"] = float(
                        re.search(
                            r"Estimated Free Energy of Binding\s*:\s*([-\d.]+)", log
                        )[1]
                    )
                scores.append(row)
    for key in {r["key"] for r in scores}:
        pool = [r for r in scores if r["key"] == key]
        for arm in ("parent", "fixed"):
            for rank, r in enumerate(
                sorted(pool, key=lambda r: (r[arm + "_kcal_mol"], r["mode"])), 1
            ):
                r[arm + "_rank"] = rank
    with (w / "frozen_pool.csv").open("w") as f:
        out = csv.DictWriter(f, scores[0], lineterminator="\n")
        out.writeheader()
        out.writerows(scores)


if __name__ == "__main__":
    main()
