"""Check K=1, nondefault redundancy radii and parallel trajectory determinism."""

import argparse
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path


def one(binary, root, out, name, cpu, runs, n, radius, common_merge):
    meta = json.loads((root / "3EML/preparation.json").read_text())
    args = [
        binary,
        "--ligand",
        root / "3UZA/ligand.pdbqt",
        "--receptor",
        root / "3EML/receptor.pdbqt",
        "--num_modes",
        str(n),
        "--min_rmsd",
        str(radius),
        "--exhaustiveness",
        str(runs),
        "--cpu",
        str(cpu),
        "--max_evals",
        "12000",
        "--seed",
        "75467",
        "--energy_range",
        "100",
        "--search_counters",
        "--out",
        out,
    ]
    if common_merge:
        args.append("--use_min_rmsd_for_merge")
    for axis, v in zip("xyz", meta["center"]):
        args += ["--center_" + axis, str(v), "--size_" + axis, "24"]
    p = subprocess.run(args, check=True, text=True, capture_output=True)
    native = []
    skip = {}
    for line in p.stdout.splitlines():
        if any(
            x in line
            for x in (
                "Search trajectory:",
                "Hunt chain:",
                "Native merge:",
                "Search counters:",
            )
        ):
            # Explicit optimizer work is allowed to change, input capacity is fixed.
            native.append(line.lstrip("*").split(" explicit_calls=")[0])
        if "Deferred hunt:" in line:
            skip = dict(re.findall(r"(\w+)=([^ ]+)", line))
    return {
        "name": name,
        "cpu": cpu,
        "runs": runs,
        "n": n,
        "radius": radius,
        "common_merge": common_merge,
        "counters": native,
        "skip": skip,
        "pool_sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "parent", "candidate", "output"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for n, radius, common_merge, runs, cpu in [
            (1, 1.0, False, 1, 1),
            (9, 3.5, False, 1, 1),
            (9, 0.0, True, 1, 1),
            (9, 1.0, False, 3, 1),
            (9, 1.0, False, 3, 2),
        ]:
            pair = [
                one(
                    getattr(a, arm),
                    a.work,
                    Path(tmp) / "pose.pdbqt",
                    arm,
                    cpu,
                    runs,
                    n,
                    radius,
                    common_merge,
                )
                for arm in ("parent", "candidate")
            ]
            if pair[0]["counters"] != pair[1]["counters"]:
                raise ValueError("MC counters or number of explicit inputs changed")
            c = pair[1]
            if n == 1 and (
                c["skip"]["replaced"] != "0"
                or c["pool_sha256"] != pair[0]["pool_sha256"]
            ):
                raise ValueError("Changed sole pose")
            if cpu == 2:
                for old, new in zip(rows[-2:], pair):
                    if (old["pool_sha256"], old["skip"]) != (
                        new["pool_sha256"],
                        new["skip"],
                    ):
                        raise ValueError("Non-deterministic parallel merge")
            rows.extend(pair)
    a.output.write_text(json.dumps(rows, indent=2) + "\n")
