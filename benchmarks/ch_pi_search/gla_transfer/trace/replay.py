"""Log the selected checkpoints without changing the frozen GLA search settings."""

import argparse
import hashlib
import os
import shlex
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from study import read_csv, timed_process, write_csv

PREFIXES = (
    "Search trajectory:",
    "Hunt chain:",
    "Pocket start:",
    "Search counters:",
    "Late pose selection:",
)


def counters(path):
    return [
        line
        for line in path.read_text().splitlines()
        if line.lstrip("*").startswith(PREFIXES)
    ]


def run(work, archive, binary, dest):
    if any(
        name in os.environ for name in ("FD_ARCHIVE_TRACE", "FASTDOCK_RESEARCH_CHPI")
    ):
        raise ValueError(
            "Clear both field and trace variables before controlled replay"
        )
    dest.mkdir(parents=True, exist_ok=True)
    rows = []
    for row in read_csv(HERE.parent / "results/costs.csv"):
        tag = "_".join(
            row[k] for k in ("panel", "query", "receptor", "seed", "cap", "search")
        )
        args = shlex.split(row["args"])
        if args[args.index("--cpu") + 1] != "1":
            raise ValueError("This observer requires one worker")
        pose = dest / (tag + ".pdbqt")
        log = dest / (tag + ".log")
        trace = dest / (tag + ".trace")
        cmd = ["env", "FD_ARCHIVE_TRACE=" + str(trace)]
        if row["search"] == "field":
            field = archive / "typing" / f"GLA_{row['query']}_{row['receptor']}.txt"
            cmd.append("FASTDOCK_RESEARCH_CHPI=" + str(field))
        args[args.index("--out") + 1] = str(pose)
        timing = timed_process([*cmd, str(binary), *args], work, log)
        original = archive / "pools" / (tag + ".pdbqt")
        if (
            pose.read_bytes() != original.read_bytes()
            or hashlib.sha256(pose.read_bytes()).hexdigest() != row["sha256"]
        ):
            raise ValueError(f"Changed released poses: {tag}")
        observations = counters(log)
        if len(observations) != 14 or observations != counters(
            archive / "pools" / (tag + ".log")
        ):
            raise ValueError(
                f"Changed counters, starting centroids, or stage counts: {tag}"
            )
        rows.append(
            dict(
                tag=tag,
                pool_sha256=row["sha256"],
                trace_sha256=hashlib.sha256(trace.read_bytes()).hexdigest(),
                matched_counter_lines=len(observations),
                baseline_wall_s=row["wall_s"],
                baseline_cpu_s=row["cpu_s"],
                **timing,
            )
        )
    write_csv(dest / "replay.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "archive", "binary", "dest"):
        p.add_argument(name, type=lambda v: Path(v).resolve())
    a = p.parse_args()
    run(a.work, a.archive, a.binary, a.dest)
