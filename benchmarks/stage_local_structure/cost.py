"""Freeze attainable four-arm CPU allocations; audit coordinates and measured comparison denominators."""

import argparse
import json
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "search_trial_panel"))
import panel


def freeze(work):
    """Use the historical CPU-only affine model; exclude rather than clip unattainable targets."""
    panel.freeze(work)
    path = work / "caps.json"
    frozen = json.loads(path.read_text())
    probes = {
        arm: panel.read(work / "pilot" / arm / "cells.csv") for arm in panel.engine.ARMS
    }
    for key, row in frozen["directions"].items():
        if "caps" not in row:
            continue
        target = row["target_cpu_s"]
        ranges = {}
        for arm in panel.engine.ARMS:
            samples = {
                int(r["cap"]): float(r["cpu_s"])
                for r in probes[arm]
                if r["query"] + "-" + r["receiver"] == key
            }
            lo, hi = samples[3000], samples[12000]
            ranges[arm] = [lo, lo + (240000 - 3000) * (hi - lo) / (12000 - 3000)]
        row["predicted_cpu_ranges_s"] = ranges
        outside = [arm for arm, (lo, hi) in ranges.items() if not lo <= target <= hi]
        if outside:
            del row["caps"]
            row["excluded"] = {"unreachable_target": outside}
    panel.dump(path, frozen)


def collect(work, output):
    caps = json.loads((work / "caps.json").read_text())["directions"]
    if any("caps" in row for row in caps.values()):
        panel.collect(work, output)
        return
    # No matched jobs were eligible. Retain zero counts and still evaluate all fixed jobs.
    panel.archive(work / "fixed", output / "fixed", 65)
    panel.archive_cost(work, output / "pilot", "pilot")
    matched = output / "matched"
    matched.mkdir()
    for name in ("cells.csv", "ranks.csv"):
        (matched / name).write_text(
            (output / "pilot" / name).read_text().splitlines()[0] + "\n"
        )
    shutil.copyfile(work / "caps.json", output / "caps.json")
    aggregate(output)


def aggregate(output):
    """Generalize longitudinal aggregation only for zero prospectively matchable directions."""
    caps = json.loads((output / "caps.json").read_text())["directions"]
    if any("caps" in row for row in caps.values()):
        panel.aggregate(output)
        return
    from evaluate import summaries
    from statistics import median, mean

    upstream = panel.read(panel.UP / "results/cells.csv")
    retry = {r["key"]: r for r in panel.read(panel.UP / "default_completion/cells.csv")}
    selected = []
    for r in upstream:
        if r["arm"] == "stable":
            continue
        r = retry.get(r["key"], r).copy()
        if r["arm"] == "candidate":
            r["arm"] = "stable"
            r["key"] = r["key"].removesuffix("-candidate") + "-stable"
        selected.append(r)
    panel.engine.csvfile(
        output / "cells.csv", selected + panel.read(output / "fixed/cells.csv")
    )
    summaries(output)
    actual = []
    totals = []
    for baseline in panel.engine.ARMS[:-1]:
        rows = [
            r
            for r in panel.read(output / "paired.csv")
            if r["cohort"] == "cost" and r["baseline"] == baseline
        ]
        for r in rows:
            r["within_25pct"] = (
                bool(r["delta_first"]) and 0.8 <= float(r["cpu_ratio"]) <= 1.25
            )
        actual.extend(rows)
        changes = [float(r["delta_first"]) for r in rows if r["within_25pct"]]
        ratios = [float(r["cpu_ratio"]) for r in rows]
        totals.append(
            dict(
                cohort="cost",
                baseline=baseline,
                scheduled=len(rows),
                within_25pct=len(changes),
                ratio_min=min(ratios),
                ratio_median=median(ratios),
                ratio_max=max(ratios),
                mean_delta_within=mean(changes) if changes else None,
                improved_within=sum(x < -1e-8 for x in changes),
                tied_within=sum(abs(x) <= 1e-8 for x in changes),
                worse_within=sum(x > 1e-8 for x in changes),
            )
        )
        totals.append(
            dict(
                cohort="matched",
                baseline=baseline,
                scheduled=0,
                within_25pct=0,
                ratio_min=None,
                ratio_median=None,
                ratio_max=None,
                mean_delta_within=None,
                improved_within=0,
                tied_within=0,
                worse_within=0,
            )
        )
    panel.engine.csvfile(output / "actual_cost_pairs.csv", actual)
    panel.engine.csvfile(output / "actual_cost_summary.csv", totals)


def check(output):
    caps = json.loads((output / "caps.json").read_text())["directions"]
    n = sum("caps" in row for row in caps.values())
    expected = {
        arm: rec["binary_sha256"]
        for arm, rec in json.loads((output / "builds.json").read_text())[
            "versions"
        ].items()
    }
    assert expected["candidate"] == (output / "fixed/binary.sha256").read_text().strip()
    for stage, count in (("fixed", 65), ("pilot", 80), ("matched", 4 * n)):
        if count:
            print(stage, panel.verify(output / stage, count))
            panel.verify_receipts(output, stage, caps)
            if stage != "fixed":
                assert expected == json.loads(
                    (output / stage / "binaries.json").read_text()
                )
        else:
            assert not panel.read(output / stage / "cells.csv")
    with tempfile.TemporaryDirectory() as temp:
        work = Path(temp)
        for arm in panel.engine.ARMS:
            out = work / "pilot" / arm
            out.mkdir(parents=True)
            panel.engine.csvfile(
                out / "cells.csv",
                [r for r in panel.read(output / "pilot/cells.csv") if r["arm"] == arm],
            )
        freeze(work)
        assert (work / "caps.json").read_bytes() == (output / "caps.json").read_bytes()
        for name in ("fixed", "matched"):
            (work / name).symlink_to(
                (output / name).resolve(), target_is_directory=True
            )
        aggregate(work)
        for name in (
            "cells.csv",
            "summary.csv",
            "paired.csv",
            "paired_summary.csv",
            "actual_cost_pairs.csv",
            "actual_cost_summary.csv",
        ):
            assert (work / name).read_bytes() == (output / name).read_bytes(), name


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    cli = parser.add_subparsers(dest="action", required=True)
    for action in ("freeze", "collect", "check"):
        cmd = cli.add_parser(action)
        if action != "check":
            cmd.add_argument("work", type=Path)
        if action != "freeze":
            cmd.add_argument("output", type=Path)
    options = vars(parser.parse_args())
    {"freeze": freeze, "collect": collect, "check": check}[options.pop("action")](
        **options
    )
