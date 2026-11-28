"""Prepare the preregistered deposits with the historical cognate pipeline."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
from urllib.error import URLError

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "geometry_metric"))
import panel

CASES = [("4OCZ", "2RU", "EPHX2"), ("6AUM", "BXV", "EPHX2"), ("4R6E", "3JD", "PARP1")]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cpu():
    return sum(r.ru_utime + r.ru_stime for r in (
        resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    work = parser.parse_args().work.resolve()
    if work.exists() and any(work.iterdir()):
        parser.error("use an empty directory to preserve the recorded inputs and costs")
    work.mkdir(parents=True, exist_ok=True)
    records = []
    for pdb, lig, target in CASES:
        start, cpu_start = time.monotonic(), cpu()
        row = dict(pdb=pdb, ligand=lig, target=target)
        try:
            panel.prepare(work, pdb, lig, target)
            row["status"] = "ready"
        except (RuntimeError, ValueError, KeyError, subprocess.CalledProcessError, URLError) as exc:
            row.update(status="excluded", reason=str(exc).replace(str(work), "WORK"))
        row.update(wall_s=time.monotonic() - start, cpu_s=cpu() - cpu_start)
        row["sha256"] = {p.name: sha(p) for p in sorted((work / pdb).iterdir()) if p.is_file()}
        records.append(row)
        print(pdb, row["status"], row.get("reason", ""), flush=True)
    (work / "manifest.json").write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    main()
