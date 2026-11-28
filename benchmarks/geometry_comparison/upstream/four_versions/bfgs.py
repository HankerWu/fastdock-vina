"""Attribute the equal-cap retrospective regression to two actual optimizer revisions."""

import argparse
import json
from pathlib import Path
import re
import statistics
import tarfile
from concurrent.futures import ProcessPoolExecutor

from evaluate import collect, read, verify
from run_panel import csvfile, digest, dock, panel

REVISIONS = ("e7031e9", "2326430")


def execute(revision, binary, out, jobs):
    raw = out / "raw"
    raw.mkdir(parents=True, exist_ok=False)
    work = [
        (revision, binary.resolve(), raw, c)
        for c in panel(out / "fixtures")
        if c[0] in ("cognate", "reciprocal")
    ]
    if len(work) != 40:
        raise ValueError("incomplete input panel")
    with ProcessPoolExecutor(jobs) as pool:
        result = list(pool.map(dock, work))
    csvfile(out / "cells.csv", [r for r, _ in result])
    csvfile(out / "ranks.csv", [m for _, records in result for m in records])
    (out / "binary.sha256").write_text(digest(binary) + "\n")


def analyze(work, core, out, core_archived):
    collect(work, out, REVISIONS)
    audit = verify(out, 80)
    datasets = [(out, REVISIONS)]
    if core_archived:
        datasets.append((core, ("parent", "stable")))
    else:
        datasets.extend((core / a, (a,)) for a in ("parent", "stable"))
    rows = []
    for path, versions in datasets:
        archive = (
            tarfile.open(path / "raw.tar.xz")
            if (path / "raw.tar.xz").exists()
            else None
        )
        try:
            for r in read(path / "cells.csv"):
                if (
                    r["cohort"] not in ("cognate", "reciprocal")
                    or r["arm"] not in versions
                ):
                    continue
                name = r["key"] + ".log"
                log = (
                    archive.extractfile(name).read().decode()
                    if archive
                    else (path / "raw" / name).read_text()
                )
                hunt = re.findall(r"\bhunt_calls=(\d+)", log)
                r["hunt_calls"] = sum(map(int, hunt)) if hunt else None
                counts = ("hunt_evals", "second_evals", "explicit_evals")
                r["total_evals"] = (
                    sum(int(r[k]) for k in counts)
                    if all(r[k] for k in counts)
                    else None
                )
                rows.append(r)
        finally:
            if archive:
                archive.close()
    key = lambda r: (r["cohort"], r["target"], r["query"], r["receiver"], r["seed"])
    index = {(key(r), r["arm"]): r for r in rows}
    pairs = []
    for cell in sorted({key(r) for r in rows}):
        for a, b in [
            ("parent", REVISIONS[0]),
            REVISIONS,
            ("2326430", "stable"),
            ("parent", "stable"),
        ]:
            old, new = index[cell, a], index[cell, b]
            result = dict(zip(("cohort", "target", "query", "receiver", "seed"), cell))
            result.update(
                previous=a, next=b, old_status=old["status"], new_status=new["status"]
            )
            for field in (
                "first",
                "best5",
                "best",
                "regret",
                "hunt_calls",
                "total_evals",
                "cpu_s",
            ):
                if (
                    old["status"] == new["status"] == "ok"
                    and old[field] is not None
                    and new[field] is not None
                ):
                    result["delta_" + field] = float(new[field]) - float(old[field])
                else:
                    result["delta_" + field] = None
            pairs.append(result)
    csvfile(out / "all_cells.csv", sorted(rows, key=lambda r: (key(r), r["arm"])))
    csvfile(out / "paired.csv", pairs)
    summaries = []
    for cohort, a, b in sorted(
        {(p["cohort"], p["previous"], p["next"]) for p in pairs}
    ):
        g = [
            p
            for p in pairs
            if (p["cohort"], p["previous"], p["next"]) == (cohort, a, b)
        ]
        vs = [p for p in g if p["delta_first"] is not None]
        r = dict(
            cohort=cohort,
            previous=a,
            next=b,
            scheduled=len(g),
            paired=len(vs),
            improved=sum(p["delta_first"] < -1e-8 for p in vs),
            worse=sum(p["delta_first"] > 1e-8 for p in vs),
            ties=sum(abs(p["delta_first"]) <= 1e-8 for p in vs),
        )
        for k in (
            "first",
            "best5",
            "best",
            "regret",
            "hunt_calls",
            "total_evals",
            "cpu_s",
        ):
            vals = [p["delta_" + k] for p in vs if p["delta_" + k] is not None]
            r["mean_delta_" + k] = statistics.mean(vals) if vals else None
        summaries.append(r)
    csvfile(out / "summary.csv", summaries)
    (out / "verification.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    subs = p.add_subparsers(dest="stage", required=True)
    e = subs.add_parser("run")
    e.add_argument("revision", choices=REVISIONS)
    e.add_argument("binary", type=Path)
    e.add_argument("out", type=Path)
    e.add_argument("--jobs", type=int, default=1)
    a = subs.add_parser("analyze")
    a.add_argument("work", type=Path)
    a.add_argument("core", type=Path)
    a.add_argument("out", type=Path)
    a.add_argument("--core-archived", action="store_true")
    args = vars(p.parse_args())
    stage = args.pop("stage")
    (execute if stage == "run" else analyze)(**args)
