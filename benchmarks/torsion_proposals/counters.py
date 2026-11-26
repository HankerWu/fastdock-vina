"""Verify native odd trajectories and measure actual optimization-call counts from logs."""

import argparse
import csv
import re
from pathlib import Path


def read(path):
    lines = path.read_text().splitlines()
    trajectories = [
        {k: int(v) for k, v in re.findall(r"(\w+)=(\d+)", line)}
        for line in lines
        if line.startswith("Search trajectory:")
    ]
    if len(trajectories) != 4:
        raise ValueError(path)
    late = next(line for line in lines if line.startswith("Search counters:"))
    return trajectories, int(re.search(r"explicit_evals=(\d+)", late)[1])


def main(work):
    rows = []
    for run in sorted(work.glob("runs_*/*_prepared.log")):
        cand, explicit = read(run)
        baseline, base_explicit = read(
            run.with_name(run.name.replace("_prepared.log", "_native.log"))
        )
        odd = all(a == b for a, b in zip(baseline[::2], cand[::2]))
        even = sum(a != b for a, b in zip(baseline[1::2], cand[1::2]))
        a, b, seed, arm = run.stem.split("_")

        def count(t):
            return sum(x["hunt_evals"] + x["second_evals"] for x in t)

        rows.append(
            dict(
                cap=run.parent.name.split("_")[1],
                ligand=a,
                receiver=b,
                seed=seed,
                odd_counters_identical=odd,
                even_counters_changed=even,
                base_search=count(baseline),
                prepared_search=count(cand),
                base_explicit=base_explicit,
                prepared_explicit=explicit,
            )
        )
    if not all(r["odd_counters_identical"] for r in rows):
        raise ValueError("native odd trajectories changed")
    with (work / "counter_audit.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    main(p.parse_args().work)
