"""Replay the original 27 protonated cells at the two fixed cost-control caps."""

import argparse
import csv
import re
import shlex
from pathlib import Path
from experiment import timed_process, write_csv


def budget(work, binary):
    with (work / "compute.csv").open() as f:
        refs = [
            r
            for r in csv.DictReader(f)
            if r["state"] == "protonated" and r["arm"] == "pro"
        ]
    if len(refs) != 27:
        raise ValueError("Incomplete baseline manifest")
    blocks = {cap: [] for cap in (12100, 24000)}
    for i, r in enumerate(refs):
        for cap in (12100, 24000) if i % 2 == 0 else (24000, 12100):
            d = work / f"budget_{cap}"
            (d / "runs").mkdir(exist_ok=True, parents=True)
            args = shlex.split(r["args"])
            args[args.index("--max_evals") + 1] = str(cap)
            pos = args.index("--out") + 1
            args[pos] = str(Path(d.name) / args[pos])
            log = (work / args[pos]).with_suffix(".log")
            times = timed_process([str(binary), *args], work, log)
            trace = log.read_text()
            searches = re.findall(r"hunt_evals=(\d+).*?second_evals=(\d+)", trace)
            explicit = re.findall(r"explicit_evals=(\d+)", trace)
            row = dict(r, **times)
            row.update(
                args=" ".join(args),
                search_evals=sum(int(a) + int(b) for a, b in searches),
                explicit_evals=sum(map(int, explicit)),
                status="ok" if (work / args[pos]).is_file() else "missing",
            )
            blocks[cap].append(row)
            write_csv(d / "compute.csv", blocks[cap])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    a = p.parse_args()
    budget(a.work.resolve(), a.binary.resolve())
