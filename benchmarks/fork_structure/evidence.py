"""Archive actual arms and verify their coordinates, commands, repeatability and costs."""

import argparse
import json
import shutil
import tarfile
from statistics import mean, median
from collections import defaultdict
from pathlib import Path

import timing
import holdout

panel = timing.panel
engine = timing.engine
evaluate = holdout.external.evaluate


def collect_external(fixtures, work, output, arms):
    evaluate.collect(work, output, arms)
    commands, exclusions = {}, []
    for arm in arms:
        commands.update(json.loads((work / arm / "commands.json").read_text()))
        exclusions.extend(
            dict(arm=arm, **r)
            for r in json.loads((work / arm / "exclusions.json").read_text())
        )
    panel.dump(output / "commands.json", commands)
    panel.dump(output / "exclusions.json", exclusions)
    panel.dump(output / "verification.json", holdout.external.verify(fixtures, output))
    if arms == engine.ARMS:
        evaluate.summaries(output)


def check_cost(path, stage):
    rows = evaluate.read(path / "cells.csv")
    with tarfile.open(path / "raw.tar.xz") as archive:
        for r in rows:
            args = json.load(archive.extractfile(r["key"] + ".argv.json"))

            def arg(k):
                return args[args.index("--" + k) + 1]

            assert (
                int(arg("seed"))
                == int(r["seed"])
                == (93529 if stage == "matched" else 74317)
            )
            assert int(arg("max_evals")) == int(r["cap"])
            assert (
                Path(arg("ligand")).parent.name == r["query"]
                and Path(arg("receptor")).parent.name == r["receiver"]
            )
            assert (
                arg("cpu") == "1"
                and arg("exhaustiveness") == "3"
                and arg("num_modes") == "9"
            )
            assert ("--search_counters" in args) == (
                r["arm"] in ("stable", "candidate")
            )
            if stage == "probes":
                assert (
                    r["cohort"] in [f"probe{i}_{r['cap']}" for i in (1, 2)]
                    and int(r["cap"]) in timing.CAPS
                )
    return evaluate.verify(path, len(rows))


def cpu_repeats(rows):
    pairs = defaultdict(list)
    for row in rows:
        pairs[row["query"], row["receiver"], row["arm"], int(row["cap"])].append(row)
    result = []
    for (query, receiver, arm, cap), reps in sorted(pairs.items()):
        assert len(reps) == 2
        rates = [float(r["cpu_s"]) for r in reps]
        ratio = max(rates) / min(rates)
        result.append(
            dict(
                query=query,
                receiver=receiver,
                arm=arm,
                cap=cap,
                cpu_1_s=rates[0],
                cpu_2_s=rates[1],
                max_min_ratio=ratio,
                repeatable=ratio <= 1.25,
                identical_pool=reps[0]["pool_sha256"] == reps[1]["pool_sha256"],
            )
        )
    return result


def collect_timing(work, out, arms):
    out.mkdir(exist_ok=True, parents=True)
    for stage in ("probes", "validation", "matched"):
        if not (work / stage).exists():
            continue
        if not any((work / stage / arm / "cells.csv").exists() for arm in arms):
            continue
        dest = out / stage
        evaluate.collect(work / stage, dest, arms)
        panel.dump(dest / "verification.json", check_cost(dest, stage))
    rows = evaluate.read(out / "probes/cells.csv")
    assert len(rows) == 60 * len(arms)
    engine.csvfile(out / "repeats.csv", cpu_repeats(rows))
    if arms == engine.ARMS:
        for name in ("caps.json", "eligible.json"):
            shutil.copyfile(work / name, out / name)
        verify_allocations(out)
        if (out / "matched/cells.csv").is_file():
            evaluate.summaries(out / "matched")
        report_actual(out / "matched", "matched")


def verify_allocations(out):
    rows = evaluate.read(out / "probes/cells.csv")
    limits = timing.freeze(rows)
    assert limits == json.loads((out / "caps.json").read_text())
    validation = (
        evaluate.read(out / "validation/cells.csv")
        if (out / "validation/cells.csv").exists()
        else []
    )
    assert len(validation) == 4 * sum("caps" in v for v in limits.values())
    eligible = timing.validated(limits, validation)
    assert eligible == json.loads((out / "eligible.json").read_text())
    matched = (
        evaluate.read(out / "matched/cells.csv")
        if (out / "matched/cells.csv").exists()
        else []
    )
    assert len(matched) == 4 * sum("caps" in v for v in eligible.values())
    for r in matched:
        assert (
            int(r["cap"])
            == eligible[r["query"] + "-" + r["receiver"]]["caps"][r["arm"]]
        )


def fixed(work, output):
    panel.archive(work, output / "fixed", 65)
    panel.verify_receipts(output, "fixed", {})
    old = evaluate.read(HERE.parent / "stage_local_structure/fixed/cells.csv")
    retained = [
        r
        for r in old
        if r["arm"] in engine.ARMS[:-1]
        and r["cohort"] in ("reciprocal", "cognate", "default", "cost")
    ]
    assert len(retained) == 195
    engine.csvfile(
        output / "cells.csv", retained + evaluate.read(output / "fixed/cells.csv")
    )
    evaluate.summaries(output)
    report_actual(output, "cost")


def report_actual(output, cohort):
    rows = (
        evaluate.read(output / "paired.csv")
        if (output / "paired.csv").is_file()
        else []
    )
    pairs = [r for r in rows if r["cohort"] == cohort]
    for r in pairs:
        r["within_25pct"] = (
            bool(r["delta_first"]) and 0.8 <= float(r["cpu_ratio"]) <= 1.25
        )
    if pairs:
        engine.csvfile(output / "actual_cost_pairs.csv", pairs)
    totals = []
    for arm in engine.ARMS[:-1]:
        selected = [r for r in pairs if r["baseline"] == arm]
        matched = [r for r in selected if r["within_25pct"]]
        deltas = [float(r["delta_first"]) for r in matched]
        ratios = [float(r["cpu_ratio"]) for r in selected if r["cpu_ratio"]]
        totals.append(
            dict(
                cohort=cohort,
                baseline=arm,
                scheduled=len(selected),
                within_25pct=len(matched),
                ratio_min=min(ratios, default=None),
                ratio_median=median(ratios) if ratios else None,
                ratio_max=max(ratios, default=None),
                mean_delta_within=mean(deltas) if deltas else None,
                improved_within=sum(d < -1e-8 for d in deltas),
                tied_within=sum(abs(d) <= 1e-8 for d in deltas),
                worse_within=sum(d > 1e-8 for d in deltas),
            )
        )
    output.mkdir(exist_ok=True)
    engine.csvfile(output / "actual_cost_summary.csv", totals)


HERE = Path(__file__).resolve().parent
if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="action", required=True)
    for action in ("external", "timing", "fixed", "check"):
        cli = sub.add_parser(action)
        if action == "external":
            cli.add_argument("fixtures", type=Path)
        if action != "check":
            cli.add_argument("work", type=Path)
        cli.add_argument("output", type=Path)
        if action in ("external", "timing"):
            cli.add_argument("--controls-only", action="store_true")
    a = p.parse_args()
    arms = engine.ARMS[:-1] if getattr(a, "controls_only", False) else engine.ARMS
    if a.action == "external":
        collect_external(a.fixtures, a.work, a.output, arms)
    elif a.action == "timing":
        collect_timing(a.work, a.output, arms)
    elif a.action == "fixed":
        fixed(a.work, a.output)
    else:
        for stage in ("probes", "validation", "matched"):
            path = a.output / stage
            if (path / "cells.csv").is_file():
                assert check_cost(path, stage) == json.loads(
                    (path / "verification.json").read_text()
                )
        if (a.output / "caps.json").exists():
            verify_allocations(a.output)
