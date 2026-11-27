"""Time complete commands in independent processes and replay both methods."""

import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

p = argparse.ArgumentParser()
for f in ("dev", "hold", "work", "probe", "vina"):
    p.add_argument("--" + f, type=Path, required=True)
a = p.parse_args()
a.work.mkdir(exist_ok=True, parents=True)
rows = []
for arm in ("lattice", "random"):
    cmd = [
        sys.executable,
        str(Path(__file__).resolve().with_name("run.py")),
        "--arms",
        arm,
    ]
    for k, v in vars(a).items():
        cmd += ["--" + k, str(v / arm if k == "work" else v)]
    t = time.perf_counter()
    with (a.work / (arm + ".log")).open("w") as f:
        sub = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
        _, status, u = os.wait4(sub.pid, 0)
        sub.returncode = os.waitstatus_to_exitcode(status)
    if sub.returncode:
        raise RuntimeError(arm)
    rows.append(
        {
            "arm": arm,
            "wall_s": time.perf_counter() - t,
            "cpu_s": u.ru_utime + u.ru_stime,
            "peak_rss_kib": u.ru_maxrss,
        }
    )
with (a.work / "cold.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
