"""Run the registered ESR1 transfer panel and verify its deposited-coordinate endpoints."""

import argparse
import itertools
import json
import re
import shutil
import sys
import subprocess
import tarfile
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
sys.path.insert(0, str(BENCH / "geometry_comparison/upstream/four_versions"))
import evaluate
import run_panel


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def structures(fixtures):
    manifest = json.loads((fixtures / "manifest.json").read_text())
    if [(r["pdb"], r["ligand"]) for r in manifest] != [
        ("1QKT", "EST"),
        ("1X7R", "GEN"),
    ]:
        raise ValueError("deposit differs from registration")
    cases, omitted = [], []
    for row in manifest:
        folder = fixtures / row["pdb"]
        for name, sha in row["sha256"].items():
            if run_panel.digest(folder / name) != sha:
                raise ValueError("input changed: " + str(folder / name))
        if row["status"] != "ready":
            omitted.append({"query": row["pdb"], "reason": row["reason"]})
            continue
        info = json.loads((folder / "metadata.json").read_text())
        lines = (folder / (row["pdb"] + ".pdb")).read_text().splitlines()
        dbrefs = [
            line
            for line in lines
            if line.startswith("DBREF ") and line[12:13] == info["chain"]
        ]
        if info["target"] != "ESR1" or not any(
            line[33:39].strip() == "P03372"
            and int(line[14:18]) == int(line[55:60])
            and int(line[20:24]) == int(line[62:67])
            for line in dbrefs
        ):
            raise ValueError("unverified UniProt numbering: " + row["pdb"])
        info["ca"] = {}
        for line in lines:
            if (
                line.startswith("ATOM  ")
                and line[21] == info["chain"]
                and line[12:16].strip() == "CA"
                and line[16] in " A"
            ):
                key = line[22:27].strip() + ":" + line[17:20]
                if key in info["ca"]:
                    raise ValueError("ambiguous C-alpha: " + key)
                info["ca"][key] = [float(line[j : j + 8]) for j in (30, 38, 46)]
        cases.append(info)
    return cases, omitted


def run(fixtures, work, arm, binary):
    out = work / arm
    raw = out / "raw"
    raw.mkdir(parents=True, exist_ok=False)
    cases, exclusions = structures(fixtures)
    commands = {}
    execute = run_panel.execute

    def record(cmd, log, timeout):
        commands[log.stem] = dict(
            argv=list(map(str, cmd)),
            timeout_s=timeout,
            max_address_bytes=run_panel.LIMIT,
        )
        return execute(cmd, log, timeout)

    run_panel.execute = record
    cells, ranks = [], []
    try:
        for query, receiver in itertools.product(cases, repeat=2):
            same = query["pdb"] == receiver["pdb"]
            try:
                ref = (
                    np.asarray(query["mapped"])
                    if same
                    else run_panel.replay.reciprocal.reference(query, receiver)
                )
            except ValueError as exc:
                exclusions.append(
                    dict(query=query["pdb"], receiver=receiver["pdb"], reason=str(exc))
                )
                continue
            for seed in run_panel.replay.SEEDS:
                spec = (
                    "cognate" if same else "reciprocal",
                    fixtures,
                    query,
                    receiver,
                    ref,
                    seed,
                    12000,
                )
                row, modes = run_panel.dock((arm, binary.resolve(), raw, spec))
                cells.append(row)
                ranks.extend(modes)
    finally:
        run_panel.execute = execute
    if cells:
        run_panel.csvfile(out / "cells.csv", cells)
        run_panel.csvfile(out / "ranks.csv", ranks)
    dump(out / "commands.json", commands)
    dump(out / "exclusions.json", exclusions)
    (out / "binary.sha256").write_text(run_panel.digest(binary) + "\n")


def collect(fixtures, work, output):
    evaluate.collect(work, output)
    evaluate.summaries(output)
    commands = {}
    exclusions = []
    for arm in run_panel.ARMS:
        commands.update(json.loads((work / arm / "commands.json").read_text()))
        exclusions += [
            dict(arm=arm, **r)
            for r in json.loads((work / arm / "exclusions.json").read_text())
        ]
    dump(output / "commands.json", commands)
    dump(output / "exclusions.json", exclusions)
    dump(output / "verification.json", verify(fixtures, output))


def verify(fixtures, output):
    with tempfile.TemporaryDirectory() as temp:
        tmp = Path(temp)
        (tmp / "raw").mkdir()
        with tarfile.open(output / "raw.tar.xz") as arc:
            arc.extractall(tmp / "raw", filter="data")
        for name in ("cells.csv", "ranks.csv"):
            shutil.copyfile(output / name, tmp / name)
        cells = evaluate.read(output / "cells.csv")
        commands = json.loads((output / "commands.json").read_text())
        if set(commands) != {row["key"] for row in cells}:
            raise ValueError("missing external command")
        for row in cells:
            key = row["key"]
            argv = commands[key]["argv"]
            positional = [v for v in argv[1:] if v != "--search_counters"]
            flags = dict(zip(positional[::2], positional[1::2], strict=True))
            if any(
                flags.get("--" + k) != v
                for k, v in {
                    "cpu": "1",
                    "seed": row["seed"],
                    "exhaustiveness": "3",
                    "max_evals": "12000",
                    "num_modes": "9",
                    "energy_range": "100",
                    "spacing": "0.375",
                }.items()
            ):
                raise ValueError("external search settings mismatch: " + key)
            if (
                Path(flags["--ligand"]).parent.name != row["query"]
                or Path(flags["--receptor"]).parent.name != row["receiver"]
            ):
                raise ValueError("external input identity mismatch: " + key)
            counters = row["arm"] in ("candidate", "stable")
            if ("--search_counters" in argv) != counters:
                raise ValueError("external counters flag mismatch: " + key)
            if counters and row["status"] == "ok":
                log = (tmp / "raw" / (key + ".log")).read_text()
                for field in (
                    "hunt_evals",
                    "second_evals",
                    "explicit_calls",
                    "explicit_evals",
                    "late_prescore",
                ):
                    values = re.findall(r"\b" + field + r"=(\d+)", log)
                    if not values or sum(map(int, values)) != int(row[field]):
                        raise ValueError("counter mismatch: " + key + " " + field)
        program = BENCH / "search_trial_nov28/verify_external.py"
        return json.loads(
            subprocess.check_output([sys.executable, program, fixtures, tmp], text=True)
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    cmd = parser.add_subparsers(dest="action", required=True)
    for stage in ("run", "collect", "check"):
        sub = cmd.add_parser(stage)
        sub.add_argument("fixtures", type=Path)
        if stage != "check":
            sub.add_argument("work", type=Path)
        if stage == "run":
            sub.add_argument("arm", choices=run_panel.ARMS)
            sub.add_argument("binary", type=Path)
        else:
            sub.add_argument("output", type=Path)
    args = vars(parser.parse_args())
    action = args.pop("action")
    if action == "check":
        print(json.dumps(verify(**args), indent=2))
    else:
        {"run": run, "collect": collect}[action](**args)
