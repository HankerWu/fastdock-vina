"""Run a source-locked four-arm cost comparison; aggregate with reused frozen controls."""

import argparse
import json
import shutil
import statistics
import sys
import tarfile
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
UP = HERE.parent / "geometry_comparison/upstream/four_versions"
sys.path.insert(0, str(UP))
import run_panel as engine
from completion import archive
from evaluate import read, summaries, verify


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def run_cell(job):
    role, binary, raw, spec = job
    # dock() hardcodes the historical official cap. A neutral label requests
    # the actual supplied cap while preserving the stock command line.
    label = "stock-official" if role == "official" else role
    observed = {}

    def record(cmd, log, timeout):
        observed["argv"] = cmd
        return execute(cmd, log, timeout)

    execute = engine.execute
    engine.execute = record
    try:
        row, ranks = engine.dock((label, binary, raw, spec))
    finally:
        engine.execute = execute
    row["arm"] = role
    dump(raw / (row["key"] + ".argv.json"), observed["argv"])
    return row, ranks


def run(args):
    args.work.mkdir(exist_ok=True, parents=True)
    output = args.work / args.stage / args.arm
    raw = output / "raw"
    raw.mkdir(parents=True)
    limits = {}
    if args.stage == "matched":
        limits = json.loads((args.work / "caps.json").read_text())["directions"]
    specs = []
    for cohort, fixtures, query, rx, ref, seed, _ in engine.panel(output / "fixtures"):
        if cohort != "reciprocal" or seed != 85231:
            continue
        direction = query["pdb"] + "-" + rx["pdb"]
        if args.stage == "matched" and "excluded" in limits[direction]:
            continue
        caps = (
            (3000, 12000)
            if args.stage == "pilot"
            else (limits[direction]["caps"][args.arm],)
        )
        for cap in caps:
            setting = "pilot" + str(cap) if args.stage == "pilot" else "matched"
            specs.append(
                (
                    setting,
                    fixtures,
                    query,
                    rx,
                    ref,
                    74317 if args.stage == "pilot" else 93529,
                    cap,
                )
            )
    jobs = [(args.arm, args.binary.resolve(), raw, s) for s in specs]
    if args.stage == "pilot" and len(jobs) != 20:
        raise ValueError("incomplete pilot")
    with ProcessPoolExecutor(args.jobs) as pool:
        data = list(pool.map(run_cell, jobs))
    engine.csvfile(output / "cells.csv", [r for r, _ in data])
    engine.csvfile(output / "ranks.csv", [r for _, ranks in data for r in ranks])
    (output / "binary.sha256").write_text(engine.digest(args.binary) + "\n")


def freeze(work):
    pilots = defaultdict(dict)
    for arm in engine.ARMS:
        for row in read(work / "pilot" / arm / "cells.csv"):
            case = row["query"] + "-" + row["receiver"]
            pilots[case].setdefault(arm, {})[int(row["cap"])] = row
    result = {}
    if len(pilots) != 10:
        raise ValueError("incomplete direction set")
    for case, arms in sorted(pilots.items()):
        if set(arms) != set(engine.ARMS) or any(
            set(p) != {3000, 12000} for p in arms.values()
        ):
            raise ValueError("missing pilot cells")
        failed = [
            a
            for a, probes in arms.items()
            if any(p["status"] != "ok" for p in probes.values())
        ]
        cpu = {
            a: (float(p[3000]["cpu_s"]), float(p[12000]["cpu_s"]))
            for a, p in arms.items()
        }
        flat = [a for a, (lo, hi) in cpu.items() if hi <= lo]
        if failed or flat:
            result[case] = {
                "excluded": {"failed_probes": failed, "nonpositive_slope": flat}
            }
            continue
        target = max(hi for lo, hi in cpu.values())
        caps = {
            a: min(240000, max(3000, round(3000 + 9000 * (target - lo) / (hi - lo))))
            for a, (lo, hi) in cpu.items()
        }
        result[case] = {"target_cpu_s": target, "caps": caps}
    dump(work / "caps.json", {"directions": result})


def archive_cost(work, out, stage):
    out.mkdir(parents=True, exist_ok=True)
    cells, ranks = [], []
    hashes = {}
    with tarfile.open(out / "raw.tar.xz", "w:xz") as tar:
        for arm in engine.ARMS:
            run = work / stage / arm
            cells.extend(read(run / "cells.csv"))
            ranks.extend(read(run / "ranks.csv"))
            hashes[arm] = (run / "binary.sha256").read_text().strip()
            for f in sorted((run / "raw").iterdir()):
                tar.add(f, arcname=f.name)
    engine.csvfile(out / "cells.csv", cells)
    engine.csvfile(out / "ranks.csv", ranks)
    dump(out / "binaries.json", hashes)
    dump(out / "verification.json", verify(out, len(cells)))


def aggregate(path):
    # The previous round's candidate is this round's stable version. Select
    # these immutable records rather than repeating their data archives.
    hist = read(UP / "results/cells.csv")
    supplement = {r["key"]: r for r in read(UP / "default_completion/cells.csv")}
    selected = []
    for r in hist:
        if r["arm"] == "stable":
            continue
        r = supplement.get(r["key"], r).copy()
        if r["arm"] == "candidate":
            r["arm"] = "stable"
            r["key"] = r["key"].removesuffix("-candidate") + "-stable"
        selected.append(r)
    selected.extend(read(path / "fixed/cells.csv"))
    selected.extend(read(path / "matched/cells.csv"))
    engine.csvfile(path / "cells.csv", selected)
    summaries(path)
    pairs = read(path / "paired.csv")
    actual = []
    for p in pairs:
        if p["cohort"] not in ("cost", "matched"):
            continue
        p["within_25pct"] = (
            bool(p["delta_first"]) and 0.8 <= float(p["cpu_ratio"]) <= 1.25
        )
        actual.append(p)
    engine.csvfile(path / "actual_cost_pairs.csv", actual)
    close = []
    for cohort in ("cost", "matched"):
        for arm in engine.ARMS[:-1]:
            pairs = [
                r for r in actual if r["cohort"] == cohort and r["baseline"] == arm
            ]
            reached = [p for p in pairs if p["within_25pct"]]
            changes = [float(p["delta_first"]) for p in reached]
            ratios = [float(p["cpu_ratio"]) for p in pairs]
            close.append(
                dict(
                    cohort=cohort,
                    baseline=arm,
                    scheduled=len(pairs),
                    within_25pct=len(reached),
                    ratio_min=min(ratios),
                    ratio_median=statistics.median(ratios),
                    ratio_max=max(ratios),
                    mean_delta_within=statistics.mean(changes) if changes else None,
                    improved_within=sum(v < -1.0e-8 for v in changes),
                    tied_within=sum(abs(v) <= 1.0e-8 for v in changes),
                    worse_within=sum(v > 1.0e-8 for v in changes),
                )
            )
    engine.csvfile(path / "actual_cost_summary.csv", close)


def collect(work, output):
    archive(work / "fixed", output / "fixed", 65)
    for stage in ("pilot", "matched"):
        archive_cost(work, output / stage, stage)
    shutil.copy2(work / "caps.json", output / "caps.json")
    aggregate(output)


def verify_receipts(path, stage, caps):
    rows = read(path / stage / "cells.csv")
    with tarfile.open(path / stage / "raw.tar.xz") as raw:
        for cell in rows:
            key = cell["key"]
            if stage in ("pilot", "matched"):
                argv = json.load(raw.extractfile(key + ".argv.json"))
                get = lambda name: argv[argv.index("--" + name) + 1]
                expected_cap = (
                    caps[cell["query"] + "-" + cell["receiver"]]["caps"][cell["arm"]]
                    if stage == "matched"
                    else int(cell["cohort"].removeprefix("pilot"))
                )
                if (
                    int(get("max_evals")) != expected_cap
                    or int(cell["cap"]) != expected_cap
                    or int(get("seed")) != (74317 if stage == "pilot" else 93529)
                    or get("cpu") != "1"
                    or get("exhaustiveness") != "3"
                    or Path(get("ligand")).parent.name != cell["query"]
                    or Path(get("receptor")).parent.name != cell["receiver"]
                    or ("--search_counters" in argv)
                    != (cell["arm"] in ("stable", "candidate"))
                ):
                    raise ValueError("command does not match frozen study: " + key)
            if cell["status"] == "ok" and cell["arm"] in ("stable", "candidate"):
                log = raw.extractfile(key + ".log").read().decode()
                for column in (
                    "hunt_evals",
                    "second_evals",
                    "explicit_calls",
                    "explicit_evals",
                    "late_prescore",
                ):
                    counts = re.findall(r"\b" + column + r"=(\d+)", log)
                    if not counts or sum(map(int, counts)) != int(cell[column]):
                        raise ValueError("counter mismatch: " + key + " " + column)


def check(path):
    controls = json.loads((UP / "builds.json").read_text())["versions"]
    expected = {
        r: controls["candidate" if r == "stable" else r]["binary_sha256"]
        for r in engine.ARMS
        if r != "candidate"
    }
    expected["candidate"] = (path / "fixed/binary.sha256").read_text().strip()
    builds = json.loads((path / "builds.json").read_text())["versions"]
    if {k: v["binary_sha256"] for k, v in builds.items()} != expected:
        raise ValueError("binary identity differs across raw runs and source record")
    for stage, n in (("fixed", 65), ("pilot", 80), ("matched", None)):
        if stage == "matched":
            caps = json.loads((path / "caps.json").read_text())["directions"]
            n = 4 * sum("caps" in row for row in caps.values())
        if (
            stage != "fixed"
            and json.loads((path / stage / "binaries.json").read_text()) != expected
        ):
            raise ValueError("engine identity mismatch")
        print(stage, verify(path / stage, n))
        verify_receipts(
            path, stage, json.loads((path / "caps.json").read_text())["directions"]
        )
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for stage in ("pilot", "matched", "fixed"):
            (tmp / stage).symlink_to((path / stage).resolve(), target_is_directory=True)
        aggregate(tmp)
        for f in (
            "cells.csv",
            "summary.csv",
            "paired.csv",
            "paired_summary.csv",
            "actual_cost_pairs.csv",
            "actual_cost_summary.csv",
        ):
            if (path / f).read_bytes() != (tmp / f).read_bytes():
                raise ValueError("aggregate mismatch: " + f)
        # Recreate the timing allocation from timing receipts only.
        (tmp / "pilot").unlink()
        for role in engine.ARMS:
            p = tmp / "pilot" / role
            p.mkdir(parents=True)
            engine.csvfile(
                p / "cells.csv",
                [r for r in read(path / "pilot/cells.csv") if r["arm"] == role],
            )
        freeze(tmp)
        if (tmp / "caps.json").read_bytes() != (path / "caps.json").read_bytes():
            raise ValueError("cost allocation mismatch")
    print(
        "Verified source-locked controls, frozen CPU allocation, scores and aggregate tables."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    cli = parser.add_subparsers(dest="action", required=True)
    for stage in ("pilot", "matched"):
        sub = cli.add_parser(stage)
        sub.add_argument("work", type=Path)
        sub.add_argument("arm", choices=engine.ARMS)
        sub.add_argument("binary", type=Path)
        sub.add_argument("--jobs", type=int, default=2)
        sub.set_defaults(stage=stage)
    cli.add_parser("freeze").add_argument("work", type=Path)
    sub = cli.add_parser("collect")
    sub.add_argument("work", type=Path)
    sub.add_argument("output", type=Path)
    cli.add_parser("check").add_argument("output", type=Path)
    options = parser.parse_args()
    if options.action == "freeze":
        freeze(options.work)
    elif options.action == "collect":
        collect(options.work, options.output)
    elif options.action == "check":
        check(options.output)
    else:
        run(options)
