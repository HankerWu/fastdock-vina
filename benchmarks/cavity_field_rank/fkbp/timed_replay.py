"""Measure process startup to the last written pool, one frozen-bank ranker per fresh process."""
import argparse
from pathlib import Path
import subprocess
import sys
import time

from frozen_assay import ARMS, CASES, HERE, cpu_seconds, table


def measure(work, output, exe):
    costs = []
    for arm in ARMS:
        root = output / arm
        (root/"pools").mkdir(parents=True, exist_ok=True)
        for name in ("plans", "capture", "affinity", "maps", "pools/base"):
            link = root/name
            if not link.exists(): link.symlink_to(work/name, target_is_directory=True)
        start, cpu = time.perf_counter(), cpu_seconds()
        subprocess.run([sys.executable, str(HERE / "frozen_assay.py"), "pipeline", str(root),
                        "--arm",arm,"--exe",str(exe)],check=True)
        elapsed = time.perf_counter()-start, cpu_seconds()-cpu
        for rec,n in CASES:
            stem = Path("pools")/arm/f"{rec}_{n}"
            for ext in ("pdbqt","xyz","stages","history"):
                path = stem.with_suffix("."+ext)
                assert (work/path).read_bytes()==(root/path).read_bytes()
        costs.append(dict(arm=arm, wall_s=elapsed[0],cpu_s=elapsed[1], verified_pools=len(CASES)))
    table(work/"end_to_end.csv",costs)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("work","output","exe"): p.add_argument(name,type=Path)
    a=p.parse_args()
    measure(a.work.resolve(),a.output.resolve(),a.exe.resolve())
