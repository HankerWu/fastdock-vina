"""Replay and audit low-cap traces, including the first sampled and minimized poses."""

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

from panel import cases, command, counts, read, sha, structural, write


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "result", "out"):
        p.add_argument(name, type=Path)
    for arm in ("parent", "reset", "ccd", "lock"):
        p.add_argument("--" + arm, type=Path)
    p.add_argument("--cap", default="3000", choices=("3000", "12000"))
    a = p.parse_args()
    raw = a.out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    stage_runs = []
    table = {case[0]: case for case in cases(a.inputs, "macro")}
    for r in read(a.result / "runs.csv"):
        if r["cap"] != a.cap or getattr(a, r["arm"], None) is None:
            continue
        key = r["key"]
        prefix = raw / key
        with tempfile.TemporaryDirectory() as tmp:
            replay = Path(tmp) / "pose.pdbqt"
            # Stage writers append models. Never use stale data from an earlier run.
            if list(raw.glob(key + ".*.pdbqt")):
                raise FileExistsError(key)
            log = subprocess.check_output(
                command(
                    table[r["pdb"]], r["seed"], r["cap"], getattr(a, r["arm"]), replay
                ),
                text=True,
                stderr=subprocess.STDOUT,
                env={**os.environ, "FASTDOCK_RING_TRACE": str(prefix)},
            )
            if sha(replay) != r["pool_sha256"] or any(
                int(r[k]) != v for k, v in counts(log).items()
            ):
                raise ValueError("Observer altered poses/work: " + key)
        for stage in ("raw", "proposal", "hunt", "second", "saved"):
            path = raw / (key + "." + stage + ".pdbqt")
            n = sum(s.startswith("MODEL") for s in path.read_text().splitlines())
            if n == 0 or (
                stage in ("hunt", "second") and n != int(r[stage + "_calls"])
            ):
                raise ValueError("Incomplete stage " + str(path))
            stage_runs.append(
                dict(
                    **{k: r[k] for k in ("pdb", "compound", "seed", "cap", "arm")},
                    key=key + "." + stage,
                    stage=stage,
                    modes=n,
                )
            )
    structural(a.inputs, a.out, stage_runs)
    write(a.out / "trace_runs.csv", stage_runs)


if __name__ == "__main__":
    main()
