"""Replay ten fixed-cap reciprocal cells with passive observers, then certify stock equivalence."""

import argparse
import json
from pathlib import Path
import tarfile

from run_panel import csvfile, digest, dock, panel
from evaluate import read


def run(arm, binary, out):
    (out / "raw").mkdir(parents=True, exist_ok=False)
    experiments = [
        cell
        for cell in panel(out / "fixtures")
        if cell[0] == "reciprocal" and cell[5] == 85231
    ]
    if len(experiments) != 10:
        raise ValueError("incomplete frozen subset")
    cells = [dock((arm, binary.resolve(), out / "raw", c))[0] for c in experiments]
    csvfile(out / "cells.csv", cells)
    (out / "binary.sha256").write_text(digest(binary) + "\n")


def verify(work, primary, out):
    out.mkdir(exist_ok=True, parents=True)
    stock = {r["key"]: r for r in read(primary / "cells.csv")}
    rows = []
    binaries = {}
    with tarfile.open(out / "raw.tar.xz", "w:xz") as archive:
        for arm in ("official", "parent"):
            arm_dir = work / arm
            binaries[arm] = (arm_dir / "binary.sha256").read_text().strip()
            for obs in read(arm_dir / "cells.csv"):
                key = obs["key"]
                base = stock[key]
                log = arm_dir / "raw" / (key + ".log")
                archive.add(log, arcname=log.name)
                pose = log.with_suffix(".pdbqt")
                if pose.exists():
                    archive.add(pose, arcname=pose.name)
                equal = (
                    base["status"] == obs["status"] == "ok"
                    and digest(pose) == base["pool_sha256"] == obs["pool_sha256"]
                )
                if equal and not all(
                    obs[k] for k in ("hunt_evals", "second_evals", "explicit_evals")
                ):
                    raise ValueError(f"missing passive count for {key}")
                rows.append(
                    dict(
                        key=key,
                        arm=arm,
                        stock_status=base["status"],
                        observer_status=obs["status"],
                        stock_sha256=base["pool_sha256"],
                        observer_sha256=obs["pool_sha256"],
                        certified=equal,
                        cpu_s=obs["cpu_s"],
                        wall_s=obs["wall_s"],
                        peak_rss_kib=obs["peak_rss_kib"],
                        **{
                            k: obs[k] if equal else None
                            for k in (
                                "hunt_evals",
                                "second_evals",
                                "explicit_evals",
                                "explicit_calls",
                                "late_prescore",
                            )
                        },
                    )
                )
    csvfile(out / "certified.csv", rows)
    (out / "binaries.json").write_text(json.dumps(binaries, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="action", required=True)
    r = sub.add_parser("run")
    r.add_argument("arm", choices=("official", "parent"))
    r.add_argument("binary", type=Path)
    r.add_argument("out", type=Path)
    v = sub.add_parser("verify")
    for k in ("work", "primary", "out"):
        v.add_argument(k, type=Path)
    args = vars(p.parse_args())
    action = args.pop("action")
    (run if action == "run" else verify)(**args)
