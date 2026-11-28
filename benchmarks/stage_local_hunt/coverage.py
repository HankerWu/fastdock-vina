"""Trace post-unblinding coverage on two RBP4 and all four ESR1 reciprocal cells."""

import argparse
import csv
import json
import lzma
import os
from pathlib import Path
import sys
import tarfile
import tempfile
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "stage_local_structure"))
from external import run_panel, structures

sys.path.insert(0, str(HERE.parent / "stage_local_attribution"))
from analyze_stages import native, observe_trajectory


def specs(tmp):
    fixed = HERE.parent / "stage_local_structure/fixed"
    ext = HERE.parent / "stage_local_structure/external"
    for job in run_panel.panel(tmp / "fixed"):
        cohort, _, query, _, _, seed, _ = job
        if cohort == "reciprocal" and query["target"] == "RBP4" and seed == 85237:
            yield fixed, job
    with tarfile.open(ext / "inputs.tar.xz") as archive:
        archive.extractall(tmp / "external", filter="data")
    cases, excluded = structures(tmp / "external")
    assert len(cases) == 2 and not excluded
    for q, r in (cases, cases[::-1]):
        ref = run_panel.replay.reciprocal.reference(q, r)
        for seed in (85231, 85237):
            yield ext, ("reciprocal", tmp / "external", q, r, ref, seed, 12000)


def audit(tmp, output, builds, rebuild):
    results, costs, transitions = [], [], []
    for panel, job in specs(tmp):
        _, inputs, q, r, ref, seed, _ = job
        coords = (inputs / q["pdb"] / "ligand.pdbqt").read_text().splitlines()
        serials = [
            int(l[6:11])
            for l in coords
            if l.startswith(("ATOM  ", "HETATM"))
            and l[77:].strip() not in {"HD", "HS", "H"}
        ]
        order = [serials.index(i) for i in q["serial"]]

        def error(row):
            xyz = np.fromstring(row[8], dtype=int, sep=" ").reshape(-1, 3)[order] / 1000
            return float(
                np.sqrt(((xyz[None] - ref) ** 2).sum(axis=-1).mean(axis=1)).min()
            )

        old = list(csv.DictReader((panel / "cells.csv").open()))
        for arm in ("stable", "candidate"):
            cell = next(
                c
                for c in old
                if c["cohort"] == "reciprocal"
                and c["query"] == q["pdb"]
                and c["receiver"] == r["pdb"]
                and int(c["seed"]) == seed
                and c["arm"] == arm
            )
            key = cell["key"]
            trace = tmp / (key + ".trace")
            packed = output / (key + ".xz")
            if builds:
                os.environ["FD_BFGS_TRACE"] = str(trace)
                os.environ.pop("FD_BFGS_SHADOW", None)
                binary = builds / arm / "build/vina"
                fresh, _ = run_panel.dock((arm, binary, output, job))
                if (
                    fresh["status"] != "ok"
                    or fresh["pool_sha256"] != cell["pool_sha256"]
                ):
                    raise ValueError(("observer changed prediction", key))
                packed.write_bytes(lzma.compress(trace.read_bytes(), preset=3))
                costs.append(
                    dict(
                        binary_sha256=run_panel.digest(binary),
                        trace_sha256=run_panel.digest(trace),
                        compressed_sha256=run_panel.digest(packed),
                        **fresh
                    )
                )
            else:
                item = next(
                    x
                    for x in json.loads((output / "costs.json").read_text())
                    if x["key"] == key
                )
                assert item["compressed_sha256"] == run_panel.digest(packed)
                trace.write_bytes(lzma.decompress(packed.read_bytes()))
                assert item["trace_sha256"] == run_panel.digest(trace)
                assert (
                    cell["pool_sha256"]
                    == item["pool_sha256"]
                    == run_panel.digest(output / (key + ".pdbqt"))
                )
            paths, explicit = native(trace)
            stages = {
                s: []
                for s in (
                    "hunt",
                    "accepted",
                    "saved",
                    "inserted",
                    "retained",
                    "merged",
                    "refined",
                )
            }
            n_evicted, n_near_evicted = 0, 0
            saved_by_state = {}
            for trajectory in paths:
                final, inserted, evicted, near_evicted = observe_trajectory(
                    trajectory, error
                )
                n_evicted += evicted
                n_near_evicted += near_evicted
                stages["inserted"].extend(inserted)
                stages["retained"].extend(v[1] for v in final.values())
                for step in trajectory:
                    h = error(step["H"][0])
                    stages["hunt"].append(h)
                    if step["A"]:
                        stages["accepted"].append(h)
                    for record in step["D"]:
                        saved_by_state[record[3]] = error(record)
                        stages["saved"].append(saved_by_state[record[3]])
            ranks = {
                row[3]: i
                for i, row in enumerate(explicit["F"][: int(cell["modes"])], 1)
            }
            before = {int(row[1]): saved_by_state[row[3]] for row in explicit["I"]}
            stages["merged"] = list(before.values())
            for record in explicit["O"]:
                after = error(record)
                stages["refined"].append(after)
                index = int(record[1])
                transitions.append(
                    dict(
                        key=key,
                        index=index,
                        pre_A=before[index],
                        post_A=after,
                        published_rank=ranks.get(record[3], ""),
                    )
                )
            published = run_panel.replay.reciprocal.measure(
                output / (key + ".pdbqt"), q["serial"], ref
            )
            assert abs(published[0] - float(cell["first"])) < 1e-9
            assert abs(min(published) - float(cell["best"])) < 1e-9
            row = dict(
                key=key,
                arm=arm,
                first=published[0],
                oracle=min(published),
                evicted=n_evicted,
                evicted_near=n_near_evicted,
            )
            stages["published"] = published
            for phase, values in stages.items():
                row[phase + "_n"] = len(values)
                row[phase + "_min"] = min(values) if values else ""
                row[phase + "_near"] = sum(v <= 2 for v in values)
            results.append(row)
    assert len(results) == 12
    if builds:
        (output / "costs.json").write_text(json.dumps(costs, indent=2) + "\n")
    for name, table in (("stages", results), ("refinement", transitions)):
        path = output / (name + ".csv")
        if builds or rebuild:
            run_panel.csvfile(path, table)
        else:
            assert list(csv.DictReader(path.open())) == [
                {k: str(v) for k, v in d.items()} for d in table
            ]
    print("Verified 12 published pools and all stage summaries")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    p.add_argument(
        "--builds",
        type=Path,
        help="fresh execution with stable and candidate observer builds",
    )
    p.add_argument(
        "--rebuild-tables",
        action="store_true",
        help="recalculate CSVs from archived traces",
    )
    a = p.parse_args()
    if a.builds:
        a.output.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory() as tmp:
        audit(
            Path(tmp).resolve(),
            a.output.resolve(),
            a.builds.resolve() if a.builds else None,
            a.rebuild_tables,
        )
