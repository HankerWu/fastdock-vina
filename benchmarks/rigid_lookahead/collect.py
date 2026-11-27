"""Collect the measured sweep, check returned ranks, and preserve its raw predictions."""

import argparse
import csv
import shutil
from collections import defaultdict
from itertools import product
from pathlib import Path
import tarfile


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def save(path, data):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=data[0], lineterminator="\n")
        w.writeheader()
        w.writerows(data)


def add_file(archive, path, name):
    info = archive.gettarinfo(str(path), str(name))
    info.uid = info.gid = info.mtime = 0
    info.uname = info.gname = ""
    with path.open("rb") as f:
        archive.addfile(info, f)


def collect(work, dest):
    costs, cells, summary = [], [], []
    process = {
        (r["stage"], int(r["budget"]), r["arm"]): r
        for r in rows(work / "process_cost.csv")
    }
    for stage, budget in product(("dev", "secondary"), (9, 24, 48)):
        with tarfile.open(dest / f"{stage}_{budget}.tar.xz", "w:xz", preset=3) as out:
            for arm in ("raw", "sham", "gradient"):
                run = work / f"{stage}_{budget}_{arm}"
                pools = defaultdict(list)
                for r in rows(run / "ranked.csv"):
                    pools[tuple(r[k] for k in ("receiver", "query", "seed"))].append(r)
                data = rows(run / "cells.csv")
                if len(data) != (12 if stage == "dev" else 8):
                    raise ValueError("incomplete panel")
                for cell in data:
                    pool = pools[tuple(cell[k] for k in ("receiver", "query", "seed"))]
                    scores = [float(r["after_rmsd"]) for r in pool]
                    if [int(r["rank"]) for r in pool] != list(
                        range(1, len(scores) + 1)
                    ):
                        raise ValueError("nonconsecutive ranks")
                    if [float(cell[k]) for k in ("top", "top5", "oracle")] != [
                        scores[0],
                        min(scores[:5]),
                        min(scores),
                    ]:
                        raise ValueError(
                            "geometry aggregation differs from ranked pool"
                        )
                    cells.append(dict(budget=budget, **cell))
                times = rows(run / "cost.csv")
                costs.extend(
                    dict(stage=stage, budget=budget, method=arm, **r) for r in times
                )
                casecost = [r for r in times if r["arm"] == arm]
                process_row = process[stage, budget, arm]
                summary.append(
                    dict(
                        stage=stage,
                        budget=budget,
                        arm=arm,
                        n=len(data),
                        **{
                            k: sum(float(r[k]) <= 2 for r in data)
                            for k in ("top", "top5", "oracle")
                        },
                        **{
                            k: sum(int(r[k]) for r in casecost)
                            for k in ("probes", "lookups", "refines", "gradients")
                        },
                        **{
                            k: process_row[k]
                            for k in ("wall_s", "cpu_s", "peak_rss_kib")
                        },
                    )
                )
                for f in sorted(run.rglob("*")):
                    relative = f.relative_to(run)
                    if f.is_file() and relative.parts[0] != "maps":
                        add_file(out, f, Path(arm) / relative)
    save(dest / "cells.csv", cells)
    save(dest / "cost.csv", costs)
    save(dest / "summary.csv", summary)


def collect_baselines(dev, secondary, dest):
    # The native development runner also measures receiver self-docks. Only
    # the six registered cross-docks belong to this comparison.
    with tarfile.open(dest / "baselines.tar.xz", "w:xz", preset=3) as out:
        for stage, folder in (("dev", dev / "runs"), ("secondary", secondary)):
            count = 0
            for pose in sorted(folder.glob("*.pdbqt")):
                q, r, *_ = pose.stem.split("_")
                if q == r:
                    continue
                for path in (pose, pose.with_suffix(".log")):
                    add_file(out, path, Path(stage) / path.name)
                count += 1
            expected = 48 if stage == "dev" else 32
            if count != expected:
                raise ValueError(f"expected {expected} native runs, found {count}")
    for filename in ("baseline_cells.csv", "baseline_cost.csv"):
        shutil.copyfile(secondary / filename, dest / filename)
    shutil.copyfile(dev / "binaries.json", dest / "binaries.json")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("dest", type=Path)
    p.add_argument(
        "--baselines",
        type=Path,
        nargs=2,
        metavar=("DEV", "SECONDARY_OUT"),
        help="also package original Vina and branch baseline predictions",
    )
    a = p.parse_args()
    a.dest.mkdir(parents=True, exist_ok=True)
    collect(a.work, a.dest)
    if a.baselines:
        collect_baselines(*a.baselines, a.dest)


if __name__ == "__main__":
    main()
