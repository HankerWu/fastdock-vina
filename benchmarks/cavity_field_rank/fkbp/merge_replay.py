"""Measure the exact AD4 merge between the per-trajectory archive and the final public pool."""
import argparse
import os
import subprocess
import time
from pathlib import Path

from frozen_assay import ARMS, CASES, cpu_seconds, sha, table


def replay(work, exe):
    tmp = work / "merge"
    tmp.mkdir(exist_ok=True)
    costs = []
    for rec, n in CASES:
        for arm in ("match", *ARMS):
            name = f"{rec}_{n}"
            orig = Path("pools") / arm / name
            root = "plans" if arm == "match" else "ranks"
            env = {k:v for k,v in os.environ.items() if not k.startswith("FD_")}
            env.update(FD_RIGID_INTERVENTION=str(work / root / f"{name}_{arm}.input"),
                       FD_RIGID_LOG=str(tmp / "probe.stages"), FD_RIGID_COORDS=str(tmp / "probe.xyz"),
                       FD_RIGID_HISTORY=str(tmp / "probe.history"))
            start, cpu = time.perf_counter(), cpu_seconds()
            subprocess.run([str(exe), "--config", str(orig.with_suffix(".conf")), "--out",
                            str(tmp / "probe.pdbqt"), "--search_counters"], cwd=work,
                           env=env, stdout=subprocess.DEVNULL, check=True)
            elapsed = time.perf_counter() - start, cpu_seconds()-cpu
            for ext in ("pdbqt","history","stages"):
                assert (work / orig.with_suffix("."+ext)).read_bytes() == (tmp / ("probe."+ext)).read_bytes()
            lines = (tmp / "probe.xyz").read_bytes().splitlines(keepends=True)
            merged = [l for l in lines if l.startswith(b"merged_ad4 ")]
            assert merged and b"".join(l for l in lines if not l.startswith(b"merged_ad4 ")) == (work / orig.with_suffix(".xyz")).read_bytes()
            (tmp / f"{name}_{arm}.xyz").write_bytes(lines[0] + b"".join(merged))
            costs.append(dict(receiver=rec, compound=n, arm=arm, merged=len(merged),
                              pool_sha256=sha(work / orig.with_suffix(".pdbqt")), exe_sha256=sha(exe),
                              wall_s=elapsed[0], cpu_s=elapsed[1]))
    table(work / "merge_cost.csv", costs)
    for p in tmp.glob("probe.*"): p.unlink()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("exe", type=Path)
    a = p.parse_args()
    replay(a.work.resolve(),a.exe.resolve())
