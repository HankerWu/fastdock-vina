"""Replay released Pro on previously disclosed ESR1/BRD4 receivers and audit all four arms."""

import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import subprocess
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import evidence

runner = evidence.holdout.external
BASE = HERE.parents[1]
FAMILIES = {
    "ESR1": (
        BASE / "stage_local_structure/external",
        BASE / "stage_local_structure/external/inputs.tar.xz",
    ),
    "BRD4": (
        BASE / "search_trial_nov28/results",
        BASE / "search_trial_nov28/external_inputs.tar.xz",
    ),
}


def brd4_inputs(inputs):
    spec = importlib.util.spec_from_file_location(
        "known_brd4", BASE / "search_trial_nov28/external.py"
    )
    brd4 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(brd4)
    manifest = json.loads((inputs / "manifest.json").read_text())
    assert [r["pdb"] for r in manifest] == ["4BJX", "4WIV", "5U2E"]
    return (
        [brd4.metadata(inputs, r) for r in manifest if r["status"] == "ready"],
        [
            dict(query=r["pdb"], reason=r["reason"])
            for r in manifest
            if r["status"] != "ready"
        ],
    )


def run(family, work, binary, output):
    prior, archive = FAMILIES[family]
    work.mkdir(parents=True, exist_ok=False)
    inputs = work / "inputs"
    with tarfile.open(archive) as tar:
        tar.extractall(inputs, filter="data")
    ids = {
        a: r["binary_sha256"]
        for a, r in json.loads((HERE.parent / "builds.json").read_text())[
            "versions"
        ].items()
    }
    assert runner.run_panel.digest(binary) == ids["candidate"]
    if family == "BRD4":
        runner.structures = brd4_inputs
    runner.run(inputs, work, "candidate", binary)
    rows = evidence.evaluate.read(prior / "cells.csv")
    ranks = evidence.evaluate.read(prior / "ranks.csv")
    historical_ids = json.loads(
        (
            prior / ("binaries.json" if family == "ESR1" else "execution.json")
        ).read_text()
    )
    if family == "BRD4":
        historical_ids = historical_ids["binary_sha256"]
    commands = json.loads(
        (
            prior / ("commands.json" if family == "ESR1" else "replay_commands.json")
        ).read_text()
    )
    with tarfile.open(prior / "raw.tar.xz") as raw:
        for arm in runner.run_panel.ARMS[:-1]:
            assert historical_ids[arm] == ids[arm]
            folder = work / arm
            data = folder / "raw"
            data.mkdir(parents=True)
            use = [r for r in rows if r["arm"] == arm]
            selected = {r["key"] for r in use}
            assert len(selected) == 8
            for key in selected:
                for ext in (".log", ".pdbqt"):
                    name = key + ext
                    member = name if family == "ESR1" else "raw/" + name
                    (data / name).write_bytes(raw.extractfile(member).read())
            runner.run_panel.csvfile(folder / "cells.csv", use)
            runner.run_panel.csvfile(
                folder / "ranks.csv", [r for r in ranks if r["key"] in selected]
            )
            (folder / "binary.sha256").write_text(ids[arm] + "\n")
            runner.dump(
                folder / "commands.json",
                {k: v for k, v in commands.items() if k in selected},
            )
            (folder / "exclusions.json").write_bytes(
                (work / "candidate/exclusions.json").read_bytes()
            )
    evidence.collect_external(inputs, work, output, runner.run_panel.ARMS)
    subprocess.run(
        [sys.executable, BASE / "stage_local_structure/directions.py", output],
        check=True,
    )


def check(family, output):
    prior, archive = FAMILIES[family]
    with tempfile.TemporaryDirectory() as tmp:
        inputs = Path(tmp)
        with tarfile.open(archive) as tar:
            tar.extractall(inputs, filter="data")
        assert runner.verify(inputs, output) == json.loads(
            (output / "verification.json").read_text()
        )
    controls = [
        r
        for r in evidence.evaluate.read(output / "cells.csv")
        if r["arm"] != "candidate"
    ]
    old = [
        r
        for r in evidence.evaluate.read(prior / "cells.csv")
        if r["arm"] != "candidate"
    ]
    assert sorted(controls, key=lambda r: r["key"]) == sorted(
        old, key=lambda r: r["key"]
    )
    with tempfile.TemporaryDirectory() as tmp:
        dest = Path(tmp)
        shutil.copy(output / "cells.csv", dest / "cells.csv")
        evidence.evaluate.summaries(dest)
        subprocess.run(
            [sys.executable, BASE / "stage_local_structure/directions.py", dest],
            check=True,
        )
        for name in (
            "summary.csv",
            "paired.csv",
            "paired_summary.csv",
            "by_direction.csv",
        ):
            assert (dest / name).read_bytes() == (output / name).read_bytes()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="action", required=True)
    for step in ("run", "check"):
        cmd = sub.add_parser(step)
        cmd.add_argument("family", choices=tuple(FAMILIES))
        if step == "run":
            cmd.add_argument("work", type=Path)
            cmd.add_argument("binary", type=Path)
        cmd.add_argument("output", type=Path)
    args = vars(p.parse_args())
    {"check": check, "run": run}[args.pop("action")](**args)
