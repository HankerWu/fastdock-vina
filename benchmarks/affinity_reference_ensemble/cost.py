"""Serial one-core cost check on the first S6K1 included by the measured table."""

import argparse
import csv
import json
import os
import re
from pathlib import Path
import subprocess
import tempfile
import time
from run import nodes


def run(a):
    ligand = next(
        r["cid"]
        for r in csv.DictReader((a.s6k / "panel.csv").open())
        if r["decision"] == "included"
    )
    center = json.loads((a.s6k / "prepared/3A60.json").read_text())["center"]
    rec = a.s6k / "prepared/3A60.pdbqt"
    lig = a.s6k / "prepared" / (ligand + ".pdbqt")
    actual = a.s6k / "pools" / f"3A60_{ligand}_271828_pro.pdbqt"
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        points = root / "nodes"
        nodes(points, 1024)
        for name in ["pro", "joint+controls", "pro_double", "upstream"]:
            if name == "joint+controls":
                cmd = [a.integrator, rec, lig, actual, points, "298.15"]
            else:
                cmd = [
                    a.upstream if name == "upstream" else a.pro,
                    "--receptor",
                    rec,
                    "--ligand",
                    lig,
                    "--seed",
                    "271828",
                    "--cpu",
                    "1",
                    "--exhaustiveness",
                    "4" if name == "pro_double" else "2",
                    "--max_evals",
                    "25000",
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--out",
                    root / (name + ".pdbqt"),
                ]
                if name != "upstream":
                    cmd += ["--pose_selection", "late", "--search_counters"]
                for k, v in zip("xyz", center):
                    cmd += ["--center_" + k, str(v), "--size_" + k, "26"]
            began = time.monotonic()
            with (root / "log").open("w") as log:
                p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
                # wait4 alone can retain the launcher RSS from before exec.
                # Read the kernel high-water mark of the actual target process.
                rss = 0
                proc = Path("/proc") / str(p.pid)
                while True:
                    try:
                        if (proc / "exe").resolve() == Path(cmd[0]).resolve():
                            status_text = (proc / "status").read_text()
                            high = re.search(r"^VmHWM:\s+(\d+) kB", status_text, re.M)
                            if high:
                                rss = max(rss, int(high[1]))
                    except (OSError, RuntimeError):
                        pass
                    pid, status, usage = os.wait4(p.pid, os.WNOHANG)
                    if pid:
                        break
                    time.sleep(.001)
                if not rss:
                    raise RuntimeError("Missing post-exec process memory measurement")
                p.returncode = os.waitstatus_to_exitcode(status)
                if p.returncode:
                    raise RuntimeError((root / "log").read_text())
            results.append(
                dict(
                    ligand=ligand,
                    seed=271828,
                    arm=name,
                    wall_s=time.monotonic() - began,
                    cpu_s=usage.ru_utime + usage.ru_stime,
                    peak_rss_KiB=rss,
                )
            )
    with (a.s6k / "serial_cost.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=results[0])
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ["s6k", "pro", "upstream", "integrator"]:
        p.add_argument(key, type=Path)
    run(p.parse_args())
