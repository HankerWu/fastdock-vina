"""Measure repeated CPU brackets before permitting independent-seed cost comparisons."""

import argparse
import json
import statistics
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "search_trial_panel"))
import panel

engine = panel.engine
CAPS = (3000, 12000, 120000)
ARMS = engine.ARMS


def run(stage, work, arm, binary, jobs):
    out = work / stage / arm
    raw = out / "raw"
    raw.mkdir(parents=True, exist_ok=False)
    specs = [
        s
        for s in engine.panel(out / "fixtures")
        if s[0] == "reciprocal" and s[-2] == 85231
    ]
    assert len(specs) == 10
    batches = []
    if stage == "probes":
        for i, caps in enumerate((CAPS, CAPS[::-1]), 1):
            for cap in caps:
                batches.append(
                    [(f"probe{i}_{cap}", *s[1:5], 74317, cap) for s in specs]
                )
    else:
        limits = json.loads(
            (
                work / ("caps.json" if stage == "validation" else "eligible.json")
            ).read_text()
        )
        batches.append(
            [
                (
                    stage,
                    *s[1:5],
                    74317 if stage == "validation" else 93529,
                    limits[s[2]["pdb"] + "-" + s[3]["pdb"]]["caps"][arm],
                )
                for s in specs
                if "caps" in limits[s[2]["pdb"] + "-" + s[3]["pdb"]]
            ]
        )
    data = []
    with ProcessPoolExecutor(jobs) as pool:
        for specs in batches:
            data.extend(
                pool.map(
                    panel.run_cell, [(arm, binary.resolve(), raw, s) for s in specs]
                )
            )
    if data:
        engine.csvfile(out / "cells.csv", [r for r, _ in data])
        engine.csvfile(out / "ranks.csv", [r for _, ranks in data for r in ranks])
    (out / "binary.sha256").write_text(engine.digest(binary) + "\n")


def freeze(rows):
    grouped = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for r in rows:
        grouped[r["query"] + "-" + r["receiver"]][r["arm"]][int(r["cap"])].append(r)
    assert len(grouped) == 10
    result = {}
    for direction, arms in sorted(grouped.items()):
        assert set(arms) == set(ARMS)
        medians, reasons, coverage = {}, [], {}
        for arm in ARMS:
            medians[arm] = []
            for cap in CAPS:
                reps = arms[arm][cap]
                assert len(reps) == 2 and len({r["cohort"] for r in reps}) == 2
                times = [float(r["cpu_s"]) for r in reps]
                if any(r["status"] != "ok" for r in reps):
                    reasons.append(f"{arm}/{cap}:failed")
                if max(times) / min(times) > 1.25:
                    reasons.append(f"{arm}/{cap}:unstable")
                medians[arm].append(statistics.median(times))
            low, middle, high = medians[arm]
            coverage[arm] = [low, high]
            if not low < middle < high:
                reasons.append(f"{arm}:nonmonotone")
        target = max(values[1] for values in medians.values())
        for arm, (lo, hi) in coverage.items():
            if not lo <= target <= hi:
                reasons.append(f"{arm}:unbracketed")
        row = dict(
            target_cpu_s=target,
            median_cpu_by_arm=medians,
            measured_ranges_s=coverage,
            exclusions=reasons,
        )
        if not reasons:
            caps = {}
            for arm, values in medians.items():
                i = 0 if target <= values[1] else 1
                caps[arm] = round(
                    CAPS[i]
                    + (target - values[i])
                    * (CAPS[i + 1] - CAPS[i])
                    / (values[i + 1] - values[i])
                )
            row["caps"] = caps
        result[direction] = row
    return result


def validated(limits, rows):
    actual = defaultdict(dict)
    for r in rows:
        actual[r["query"] + "-" + r["receiver"]][r["arm"]] = r
    for direction, row in limits.items():
        if "caps" not in row:
            continue
        assert set(actual[direction]) == set(ARMS)
        ratios = {
            a: float(r["cpu_s"]) / row["target_cpu_s"]
            for a, r in actual[direction].items()
        }
        row["validation_ratios"] = ratios
        rejected = [
            a
            for a, v in ratios.items()
            if not 0.8 <= v <= 1.25 or actual[direction][a]["status"] != "ok"
        ]
        if rejected:
            row["exclusions"] = [a + ":validation" for a in rejected]
            row["rejected_caps"] = row.pop("caps")
    return limits


def all_rows(work, stage):
    return [
        r
        for arm in ARMS
        for r in (
            panel.read(work / stage / arm / "cells.csv")
            if (work / stage / arm / "cells.csv").exists()
            else []
        )
    ]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="action", required=True)
    for stage in ("probes", "validation", "matched", "freeze", "eligible"):
        cmd = sub.add_parser(stage)
        cmd.add_argument("work", type=Path)
        if stage not in ("freeze", "eligible"):
            cmd.add_argument("arm", choices=ARMS)
            cmd.add_argument("binary", type=Path)
            cmd.add_argument("--jobs", type=int, default=3)
    a = p.parse_args()
    if a.action == "freeze":
        panel.dump(a.work / "caps.json", freeze(all_rows(a.work, "probes")))
    elif a.action == "eligible":
        panel.dump(
            a.work / "eligible.json",
            validated(
                json.loads((a.work / "caps.json").read_text()),
                all_rows(a.work, "validation"),
            ),
        )
    else:
        run(a.action, a.work, a.arm, a.binary, a.jobs)
