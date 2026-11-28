"""Execute the registered new BRD4 panel; retain all preparation and alignment exclusions."""

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "geometry_comparison/upstream/four_versions"))
import evaluate
import run_panel


def metadata(fixtures, rec):
    case = fixtures / rec["pdb"]
    for name, digest in rec["sha256"].items():
        if run_panel.digest(case / name) != digest:
            raise ValueError(f"changed input {case.name}/{name}")
    meta = json.loads((case / "metadata.json").read_text())
    if meta["target"] != "BRD4":
        raise ValueError("family identity mismatch")
    pdb = (case / (case.name + ".pdb")).read_text().splitlines()
    if not any(
        l.startswith("DBREF ") and l[12:13] == meta["chain"] and "O60885" in l
        for l in pdb
    ):
        raise ValueError(f"cannot confirm BRD4 UniProt mapping for {case.name}")
    ca = {}
    for l in pdb:
        if not (
            l.startswith("ATOM  ")
            and l[21] == meta["chain"]
            and l[12:16].strip() == "CA"
            and l[16] in " A"
        ):
            continue
        key = l[22:27].strip() + ":" + l[17:20]
        if key in ca:
            raise ValueError(f"ambiguous C-alpha {key}")
        ca[key] = [float(l[j : j + 8]) for j in (30, 38, 46)]
    meta["ca"] = ca
    return meta


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fixtures", type=Path)
    p.add_argument("output", type=Path)
    for arm in run_panel.ARMS:
        p.add_argument("--" + arm, required=True, type=Path)
    a = p.parse_args()
    fixtures, out = a.fixtures.resolve(), a.output.resolve()
    raw = out / "raw"
    raw.mkdir(parents=True, exist_ok=False)
    binaries = {arm: getattr(a, arm).resolve(strict=True) for arm in run_panel.ARMS}
    manifest = json.loads((fixtures / "manifest.json").read_text())
    assert [r["pdb"] for r in manifest] == ["4BJX", "4WIV", "5U2E"]
    active, exclusions = [], []
    for rec in manifest:
        if rec["status"] != "ready":
            exclusions.append(dict(query=rec["pdb"], reason=rec["reason"]))
            continue
        try:
            active.append(metadata(fixtures, rec))
        except ValueError as exc:
            exclusions.append(dict(query=rec["pdb"], reason=str(exc)))
    # Capture argv at process dispatch while retaining the historical four-arm
    # runner's timing, timeout and memory accounting.
    commands = {}
    execute = run_panel.execute

    def record(cmd, log, timeout):
        commands[log.stem] = dict(
            argv=list(map(str, cmd)),
            timeout_s=timeout,
            address_limit_bytes=run_panel.LIMIT,
        )
        return execute(cmd, log, timeout)

    run_panel.execute = record
    cells, ranks = [], []
    for query, receiver in itertools.product(active, repeat=2):
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
        cohort = "cognate" if same else "reciprocal"
        for i, seed in enumerate(run_panel.replay.SEEDS):
            job = (cohort, fixtures, query, receiver, ref, seed, 12000)
            for arm in (run_panel.ARMS if i == 0 else reversed(run_panel.ARMS)):
                cell, modes = run_panel.dock((arm, binaries[arm], raw, job))
                cells.append(cell)
                ranks.extend(modes)
    (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
    receipt = dict(
        exclusions=exclusions,
        cases=len(cells),
        binary_sha256={arm: run_panel.digest(f) for arm, f in binaries.items()},
    )
    (out / "execution.json").write_text(json.dumps(receipt, indent=2) + "\n")
    if cells:
        run_panel.csvfile(out / "cells.csv", cells)
        run_panel.csvfile(out / "ranks.csv", ranks)
        evaluate.summaries(out)


if __name__ == "__main__":
    main()
