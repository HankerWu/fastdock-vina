"""Rebuild maps in a fresh single-threaded process for each registered arm and budget."""

import argparse
import csv
import itertools
import os
from pathlib import Path
import subprocess
import sys
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("dev", "secondary", "work", "probe", "vina"):
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    work = a.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with (work / "process_cost.csv").open("w", newline="") as table:
        writer = csv.writer(table, lineterminator="\n")
        writer.writerow(["stage", "budget", "arm", "wall_s", "cpu_s", "peak_rss_kib"])
        for stage, budget, arm in itertools.product(
            ("dev", "secondary"), (9, 24, 48), ("raw", "sham", "gradient")
        ):
            out = work / f"{stage}_{budget}_{arm}"
            out.mkdir()
            cmd = [
                sys.executable,
                str(Path(__file__).with_name("panel.py")),
                "--" + stage,
                str(getattr(a, stage).resolve()),
                "--work",
                str(out),
                "--probe",
                str(a.probe.resolve()),
                "--vina",
                str(a.vina.resolve()),
                "--arms",
                arm,
                "--budget",
                str(budget),
            ]
            t = time.perf_counter()
            with (out / "process.log").open("w") as log:
                child = subprocess.Popen(
                    cmd, env=env, stdout=log, stderr=subprocess.STDOUT
                )
                _, status, usage = os.wait4(child.pid, 0)
                child.returncode = os.waitstatus_to_exitcode(status)
            if child.returncode:
                raise RuntimeError(f"failed process; see {out / 'process.log'}")
            writer.writerow(
                [
                    stage,
                    budget,
                    arm,
                    time.perf_counter() - t,
                    usage.ru_utime + usage.ru_stime,
                    usage.ru_maxrss,
                ]
            )
            table.flush()


if __name__ == "__main__":
    main()
