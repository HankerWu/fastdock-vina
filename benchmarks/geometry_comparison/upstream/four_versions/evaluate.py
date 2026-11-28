"""Verify released coordinates independently, and summarize the fixed four-version comparison."""

import argparse
import csv
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
import tarfile

import numpy as np

from run_panel import ARMS, HERE, PRO_COUNTERS, ROOT, csvfile
from verify import references

sys.path.insert(0, str(ROOT / "benchmarks/prediction_metrics"))
from pose import quantiles


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def verify(out, expected_rows=4 * 65):
    frozen = json.loads((HERE / "inputs.json").read_text())
    metadata = {}
    for mode, spec in [("reciprocal", frozen), ("cognate", frozen["cognate"])]:
        with tarfile.open(ROOT / spec["archive"]) as tar:
            for pdb in spec["inputs"]:
                name = "metadata.json" if mode == "cognate" else "preparation.json"
                member = next(
                    m
                    for m in tar.getmembers()
                    if m.name.removeprefix("./") == f"{pdb}/{name}"
                )
                data = tar.extractfile(member).read()
                if hashlib.sha256(data).hexdigest() != spec["inputs"][pdb][name]:
                    raise ValueError(f"changed experimental reference: {pdb}")
                metadata[mode, pdb] = json.loads(data)
    cells, ranks = read(out / "cells.csv"), read(out / "ranks.csv")
    if len(cells) != expected_rows or len({r["key"] for r in cells}) != len(cells):
        raise ValueError("incomplete or duplicated panel")
    modes = defaultdict(list)
    for rank in ranks:
        modes[rank["key"]].append(rank)
    run_count = 0
    maximum = 0.0
    expected = {}
    with tarfile.open(out / "raw.tar.xz") as tar:
        for row in cells:
            key = row["key"]
            log = tar.extractfile(key + ".log").read().decode()
            if row["status"] != "ok":
                if modes[key] or row["first"]:
                    raise ValueError(f"failed run scored: {key}")
                continue
            mode = "cognate" if row["cohort"] == "cognate" else "reciprocal"
            query, receiver = [metadata[mode, row[k]] for k in ("query", "receiver")]
            pair = (mode, row["query"], row["receiver"])
            if pair not in expected:
                expected[pair] = references(query, receiver)
            serial = query["serial"]
            data = tar.extractfile(key + ".pdbqt").read()
            if hashlib.sha256(data).hexdigest() != row["pool_sha256"]:
                raise ValueError(f"output hash differs: {key}")
            poses = []
            scores = []
            for model in data.decode().split("ENDMDL"):
                coords = {
                    int(s[6:11]): tuple(float(s[a : a + 8]) for a in (30, 38, 46))
                    for s in model.splitlines()
                    if s.startswith(("ATOM  ", "HETATM"))
                }
                if not coords:
                    continue
                scores.append(
                    float(
                        next(
                            s
                            for s in model.splitlines()
                            if s.startswith("REMARK VINA RESULT:")
                        ).split()[3]
                    )
                )
                observed = [coords[i] for i in serial]
                poses.append(
                    min(
                        math.sqrt(
                            math.fsum(
                                math.dist(x, y) ** 2 for x, y in zip(observed, ref)
                            )
                            / len(serial)
                        )
                        for ref in expected[pair]
                    )
                )
            stored = modes[key]
            if scores != [float(m["score"]) for m in stored] or list(
                range(1, len(scores) + 1)
            ) != [int(m["rank"]) for m in stored]:
                raise ValueError(f"rank mismatch: {key}")
            if int(row["modes"]) != len(poses) or scores[0] != float(
                row["native_score"]
            ):
                raise ValueError(f"pool mismatch: {key}")
            observed = poses + [
                poses[0],
                min(poses[:5]),
                min(poses),
                poses[0] - min(poses),
            ]
            recorded = [float(m["rmsd"]) for m in stored] + [
                float(row[k]) for k in ("first", "best5", "best", "regret")
            ]
            maximum = max(maximum, max(abs(a - b) for a, b in zip(observed, recorded)))
            if row["arm"] in PRO_COUNTERS and "Search counters: " not in log:
                raise ValueError(f"missing counters: {key}")
            run_count += 1
    if maximum > 1e-9:
        raise ValueError(f"RMSD replay mismatch: {maximum}")
    return dict(
        scheduled=len(cells),
        verified_valid=run_count,
        verified_ranks=len(ranks),
        maximum_RMSD_difference_A=maximum,
    )


def summaries(out):
    cells = read(out / "cells.csv")
    groups = defaultdict(list)
    by_case = defaultdict(dict)
    for r in cells:
        groups[r["cohort"], r["arm"], "ALL"].append(r)
        groups[r["cohort"], r["arm"], r["target"]].append(r)
        case = tuple(r[k] for k in ("cohort", "target", "query", "receiver", "seed"))
        by_case[case][r["arm"]] = r
    rows = []
    for (cohort, arm, target), group in sorted(groups.items()):
        valid = [r for r in group if r["status"] == "ok"]
        vals = np.array([float(r["first"]) for r in valid])
        entry = dict(
            cohort=cohort,
            arm=arm,
            target=target,
            scheduled=len(group),
            valid=len(valid),
            empty=sum(r["status"] == "empty" for r in group),
            failed=sum(r["status"] not in ("ok", "empty") for r in group),
            serious_gt5=sum(vals > 5),
            cpu_s=sum(float(r["cpu_s"]) for r in group),
            wall_s_sum=sum(float(r["wall_s"]) for r in group),
            peak_rss_kib=max(int(r["peak_rss_kib"]) for r in group),
        )
        probs = (0.25, 0.5, 0.75, 0.9, 0.95)
        qs = (
            quantiles(vals, np.full(len(vals), 1 / len(vals)), probs)
            if len(vals)
            else [None] * 5
        )
        entry.update(zip(("q1", "median", "q3", "p90", "p95"), qs))
        for metric in ("first", "best5", "best", "regret"):
            entry["mean_" + metric] = (
                statistics.mean(float(r[metric]) for r in valid) if valid else None
            )
        rows.append(entry)
    csvfile(out / "summary.csv", rows)
    details = []
    for case, arms in sorted(by_case.items()):
        if set(arms) != set(ARMS):
            raise ValueError(f"incomplete paired outputs: {case}")
        current = arms["candidate"]
        for name in ARMS[:-1]:
            other = arms[name]
            valid = other["status"] == current["status"] == "ok"
            r = dict(zip(("cohort", "target", "query", "receiver", "seed"), case))
            r.update(
                baseline=name,
                baseline_status=other["status"],
                candidate_status=current["status"],
                cpu_ratio=float(current["cpu_s"]) / float(other["cpu_s"])
                if float(other["cpu_s"])
                else None,
                pool_identical=(current["pool_sha256"] == other["pool_sha256"])
                if valid
                else None,
            )
            for metric in ("first", "best5", "best", "regret"):
                r["delta_" + metric] = (
                    float(current[metric]) - float(other[metric]) if valid else None
                )
            details.append(r)
    csvfile(out / "paired.csv", details)
    paired = defaultdict(list)
    for d in details:
        for scope in ("ALL", d["target"]):
            paired[d["cohort"], d["baseline"], scope].append(d)
    summary = []
    for (cohort, base, target), all_pairs in sorted(paired.items()):
        valid = [r for r in all_pairs if r["delta_first"] is not None]
        values = [r["delta_first"] for r in valid]
        ratios = [r["cpu_ratio"] for r in valid]
        summary.append(
            dict(
                cohort=cohort,
                baseline=base,
                target=target,
                n=len(all_pairs),
                paired=len(valid),
                improved=sum(v < -1e-8 for v in values),
                tied=sum(abs(v) <= 1e-8 for v in values),
                worse=sum(v > 1e-8 for v in values),
                byte_equal=sum(r["pool_identical"] for r in valid),
                mean_delta=statistics.mean(values) if values else None,
                median_delta=statistics.median(values) if values else None,
                min_delta=min(values, default=None),
                max_delta=max(values, default=None),
                min_cpu_ratio=min(ratios, default=None),
                median_cpu_ratio=statistics.median(ratios) if ratios else None,
                max_cpu_ratio=max(ratios, default=None),
            )
        )
    csvfile(out / "paired_summary.csv", summary)


def collect(work, out, arms=ARMS):
    out.mkdir(exist_ok=True, parents=True)
    with tarfile.open(out / "raw.tar.xz", "w:xz") as tar:
        for arm in arms:
            for path in sorted((work / arm / "raw").iterdir()):
                tar.add(path, arcname=path.name)
    for name in ("cells.csv", "ranks.csv"):
        rows = [r for a in arms for r in read(work / a / name)]
        csvfile(
            out / name, sorted(rows, key=lambda r: (r["key"], int(r.get("rank", 0))))
        )
    hashes = {a: (work / a / "binary.sha256").read_text().strip() for a in arms}
    (out / "binaries.json").write_text(json.dumps(hashes, indent=2) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out", type=Path)
    ap.add_argument(
        "--collect",
        type=Path,
        help="combine freshly completed per-engine work directories",
    )
    args = ap.parse_args()
    if args.collect:
        collect(args.collect, args.out)
    result = verify(args.out)
    summaries(args.out)
    (args.out / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
