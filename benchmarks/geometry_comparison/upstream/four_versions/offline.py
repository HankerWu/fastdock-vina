"""Read-only replay of every released geometry, stock counter witness, and coverage ledger."""

import hashlib
import json
from pathlib import Path
import re
import shutil
import tarfile
import tempfile

from completion import combine
from evaluate import read, summaries, verify
from run_panel import HERE


def main():
    for group, count in [
        ("results", 260),
        ("intermediate", 65),
        ("bfgs_results", 80),
        ("default_completion", 1),
    ]:
        stored = HERE / group
        audit = verify(stored, count)
        certified = json.loads((stored / "verification.json").read_text())
        if (
            any(
                audit[k] != v
                for k, v in certified.items()
                if k != "maximum_RMSD_difference_A"
            )
            or audit["maximum_RMSD_difference_A"] > 1e-9
        ):
            raise ValueError(f"independent replay failed for {group}")
        print(group, json.dumps(audit))
    primary = {r["key"]: r for r in read(HERE / "results/cells.csv")}
    counters = HERE / "counter_results"
    records = read(counters / "certified.csv")
    if len(records) != 20 or len({r["key"] for r in records}) != 20:
        raise ValueError("counter panel incomplete")
    with tarfile.open(counters / "raw.tar.xz") as tar:
        for r in records:
            key = r["key"]
            s = primary[key]
            pose = tar.extractfile(key + ".pdbqt").read()
            h = hashlib.sha256(pose).hexdigest()
            if not (
                r["certified"] == "True"
                and h == s["pool_sha256"] == r["stock_sha256"] == r["observer_sha256"]
            ):
                raise ValueError(f"counter output mismatch: {key}")
            log = tar.extractfile(key + ".log").read().decode()
            for k in (
                "hunt_evals",
                "second_evals",
                "explicit_evals",
                "explicit_calls",
                "late_prescore",
            ):
                numbers = list(map(int, re.findall(r"\b" + k + r"=(\d+)", log)))
                if not numbers or sum(numbers) != int(r[k]):
                    raise ValueError(f"counter mismatch: {key} {k}")
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp)
        shutil.copy2(HERE / "results/cells.csv", path / "cells.csv")
        summaries(path)
        for k in ("summary.csv", "paired.csv", "paired_summary.csv"):
            if read(path / k) != read(HERE / "results" / k):
                raise ValueError(f"aggregate mismatch: {k}")
        combine(HERE / "results", path, HERE / "default_completion")
        if read(path / "completed_defaults.csv") != read(
            HERE / "default_completion/completed_defaults.csv"
        ):
            raise ValueError("default completion mismatch")
    print(
        "Verified 20 stock-identical counter replays, all primary summaries and five complete official defaults."
    )


if __name__ == "__main__":
    main()
