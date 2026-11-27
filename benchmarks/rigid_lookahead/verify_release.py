"""Check the final release CLI against the original explicit-refinement observer."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


def verify(base, meta, run, binary, temp):
    with (run / "poses.csv").open() as f:
        rows = [r for r in csv.DictReader(f) if r["arm"] == "gradient"]
    count = 0
    for row in rows:
        rec, query, seed = (row[k] for k in ("receiver", "query", "seed"))
        candidate = (
            run / "runs" / f"{rec}_{query}_{seed}_gradient" / (row["slot"] + ".pdbqt")
        )
        pose = temp / "local.pdbqt"
        cmd = [
            binary,
            "--receptor",
            base / rec / "receptor.pdbqt",
            "--ligand",
            candidate,
            "--cpu",
            "1",
            "--seed",
            seed,
            "--local_only",
            "--out",
            pose,
        ]
        for axis, val in zip("xyz", meta[rec]["center"]):
            cmd.extend(["--center_" + axis, str(val), "--size_" + axis, "24"])
        subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
        if pose.read_bytes() != Path(str(candidate) + ".refined").read_bytes():
            raise ValueError(f"release refinement differs: {candidate}")
        count += 1
    return count


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in (
        "vina",
        "dev",
        "secondary",
        "sweep",
        "transfer_input",
        "transfer_run",
        "output",
    ):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    meta = json.loads((a.dev / "metadata.json").read_text())
    stages = [("dev", a.dev, meta, a.sweep / "dev_9_gradient", 108)]
    for name, base, run, n in [
        ("secondary", a.secondary, a.sweep / "secondary_9_gradient", 72),
        ("transfer", a.transfer_input, a.transfer_run, 72),
    ]:
        meta = {
            f.parent.name: json.loads(f.read_text())
            for f in base.glob("*/preparation.json")
        }
        stages.append((name, base, meta, run, n))
    counts = {}
    with tempfile.TemporaryDirectory() as w:
        for name, base, meta, run, n in stages:
            counts[name] = verify(base, meta, run, a.vina.resolve(), Path(w))
            if counts[name] != n:
                raise ValueError(f"{name}: expected {n} poses")
    a.output.write_text(
        json.dumps(
            {
                "byte_identical_to_observer": counts,
                "release_sha256": hashlib.sha256(a.vina.read_bytes()).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
