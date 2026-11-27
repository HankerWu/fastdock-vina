"""Post-hoc 10% lower-cost sensitivity, including a parent control at the same cap."""

import argparse
import itertools
import os
import random
import re
from decimal import Decimal
from pathlib import Path

from replay import BENCH, cpu, digest
from run import COUNT, cells, dock, read, write
from scipy.stats import kendalltau, spearmanr


def run(args):
    work = args.work
    (work / "raw").mkdir(parents=True, exist_ok=True)
    os.sched_setaffinity(0, {args.cpu})
    entries = cells(work, "Gal3")
    random.Random(309).shuffle(entries)
    arms = [("parent", 12000), ("candidate", 10800), ("parent", 10800)]
    data = []
    for index, cell in enumerate(entries):
        offset = index % len(arms)
        for arm, cap in arms[offset:] + arms[:offset]:
            before = cpu()
            row = dock(cell, work, getattr(args, arm), arm, 391781, cap, True)
            row["cpu_s"] = cpu() - before
            raw = (work / "raw" / (row["key"] + ".pdbqt")).read_bytes()
            row["sha256"] = digest(raw)
            if cap == 12000:
                assert (
                    raw == (args.frozen / "raw" / (row["key"] + ".pdbqt")).read_bytes()
                )
            data.append(row)
    write(work / "runs.csv", data)
    evaluate(work)


def evaluate(work):
    data = read(work / "runs.csv")
    labels = {
        r["compound"]: r
        for r in read(BENCH / "pair_gradient/galectin_external/manifest.csv")
    }
    kd = {c: Decimal(r["kd_uM"]) for c, r in labels.items()}
    ids = sorted(kd, key=kd.get)
    edges = [(a, b) for a, b in itertools.combinations(ids, 2) if kd[b] >= 3 * kd[a]]
    ranked = []
    for arm, cap in (("parent", 12000), ("parent", 10800), ("candidate", 10800)):
        subset = [r for r in data if r["arm"] == arm and int(r["cap"]) == cap]
        scores = {}
        for r in subset:
            raw = (work / "raw" / (r["key"] + ".pdbqt")).read_bytes()
            assert digest(raw) == r["sha256"]
            log = (work / "raw" / (r["key"] + ".log")).read_text()
            for k in COUNT:
                assert re.findall(r"\b" + k + r"=(\d+)", log) == [r[k]]
            assert sum(int(r[k]) for k in COUNT[:3]) == int(r["joint"])
            scores[r["compound"], r["receptor"]] = Decimal(
                raw.splitlines()[1].split()[3].decode()
            )
        for rec in ("1KJL", "3T1L", "mean"):
            recs = ["1KJL", "3T1L"] if rec == "mean" else [rec]
            sc = {c: sum(scores[c, r] for r in recs) for c in ids}
            rows = [r for r in subset if r["receptor"] in recs]

            def credit(i, j, scores=sc):
                return (scores[i] < scores[j]) + 0.5 * (scores[i] == scores[j])

            measured = [float(kd[c]) for c in ids]
            predicted = [float(sc[c]) for c in ids]
            ranked.append(
                {
                    "arm": arm,
                    "cap": cap,
                    "receptor": rec,
                    "credit": sum(credit(i, j) for i, j in edges),
                    "matched": sum(
                        credit(i, j)
                        for i, j in edges
                        if labels[i]["heavy"] == labels[j]["heavy"]
                        and labels[i]["charge"] == labels[j]["charge"]
                    ),
                    "spearman": spearmanr(predicted, measured).statistic,
                    "kendall": kendalltau(predicted, measured).statistic,
                    "joint": sum(int(r["joint"]) for r in rows),
                    "wall_s": sum(int(r["elapsed_ns"]) for r in rows) / 1e9,
                    "cpu_s": sum(float(r["cpu_s"]) for r in rows),
                }
            )
    write(work / "summary.csv", ranked)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("run", "audit"))
    p.add_argument("work", type=Path)
    for name in ("frozen", "parent", "candidate"):
        p.add_argument("--" + name, type=Path)
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    a = p.parse_args()
    if a.action == "run":
        run(a)
    else:
        evaluate(a.work)
