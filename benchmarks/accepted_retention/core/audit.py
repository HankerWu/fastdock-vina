"""Recheck stored pose digests, budget, one-slot membership and per-case endpoints."""

import argparse
import csv
import hashlib
import statistics
import tarfile
from collections import defaultdict
from pathlib import Path


def load(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def dump(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def audit(folder):
    cells = load(folder / "cells.csv")
    by_case = defaultdict(dict)
    by_arm = defaultdict(list)
    modes = defaultdict(list)
    for mode in load(folder / "modes.csv"):
        modes[mode["key"]].append(float(mode["rmsd_A"]))
    with tarfile.open(folder / "runs.tar.xz") as tar:
        for cell in cells:
            key = cell["key"]
            raw = tar.extractfile(f"runs/{key}.pdbqt").read()
            if hashlib.sha256(raw).hexdigest() != cell["pool_sha256"]:
                raise ValueError(f"Pool digest mismatch for {key}")
            xyz = modes[key]
            if [xyz[0], min(xyz[:5]), min(xyz)] != [
                float(cell[x]) for x in ("top1_A", "best5_A", "oracle_A")
            ]:
                raise ValueError(f"Endpoints disagree with per-mode rows for {key}")
            by_case[key.rsplit("_", 1)[0]][cell["arm"]] = cell
            by_arm[cell["arm"], cell["target"], int(cell["cap"])].append(cell)
    if len(by_case) != 48 or len(cells) != 144:
        raise ValueError("Incomplete discovery")
    conserved = (
        "hunt_evals",
        "hunt_calls",
        "second_evals",
        "second_calls",
        "explicit_inputs",
    )
    for case, g in by_case.items():
        base = g["parent"]
        for name in ("candidate", "ungated"):
            r = g[name]
            if (
                not all(r[x] == base[x] for x in conserved)
                or r["slots"] != base["explicit_inputs"]
            ):
                raise ValueError("Changed native search or slot budget")
            if r["replaced"] == "0" and r["pool_sha256"] != base["pool_sha256"]:
                raise ValueError("An unchanged cell changed output")
            if (
                name == "candidate"
                and r["replaced"] == "1"
                and float(r["pending_energy"]) >= float(r["previous_energy"])
            ):
                raise ValueError(
                    "A worse or equal soft-cap energy passed the strict gate"
                )
    stages = defaultdict(dict)
    for row in load(folder / "stages.csv"):
        case, arm = row["key"].rsplit("_", 1)
        stages[case, arm][row["stage"], int(row["input_index"])] = row
    stage_summary = []
    for case, g in by_case.items():
        if g["candidate"]["replaced"] != "1":
            continue
        slot = int(g["candidate"]["slot"])
        parent, trial = stages[case, "parent"], stages[case, "candidate"]
        if (
            len(parent) != 3 * int(g["parent"]["explicit_inputs"])
            or parent.keys() != trial.keys()
        ):
            raise ValueError("Incomplete passive trace")
        for (stage, i), p in parent.items():
            if (stage == "search" or slot != i) and p["coords_sha256"] != trial[
                stage, i
            ]["coords_sha256"]:
                raise ValueError(
                    "Unaffected slot changed or a change preceded native preselection"
                )
        stage_summary.append(
            {
                "case": case,
                "slot": slot,
                "displaced_in_A": parent["before", slot]["rmsd_A"],
                "displaced_out_A": parent["after", slot]["rmsd_A"],
                "incoming_in_A": trial["before", slot]["rmsd_A"],
                "incoming_out_A": trial["after", slot]["rmsd_A"],
            }
        )
    dump(folder / "changed_slots.csv", sorted(stage_summary, key=lambda x: x["case"]))
    result = []
    for (arm, target, cap), group in sorted(by_arm.items()):

        def count(column, rows=group):
            return sum(int(r[column] or 0) for r in rows)

        def mean(column, rows=group):
            return statistics.mean(float(r[column]) for r in rows)

        def hits(column, rows=group):
            return sum(float(r[column]) <= 2 for r in rows)

        result.append(
            {
                "arm": arm,
                "target": target,
                "cap": cap,
                "n": len(group),
                "replaced": count("replaced"),
                "top1_hits": hits("top1_A"),
                "top5_hits": hits("best5_A"),
                "oracle_hits": hits("oracle_A"),
                "top1_mean_A": mean("top1_A"),
                "oracle_mean_A": mean("oracle_A"),
                "top1_clash_pairs": count("top1_clash_pairs"),
                "calls": sum(
                    count(x) for x in ("hunt_evals", "second_evals", "explicit_evals")
                ),
                "seconds": sum(float(r["seconds"]) for r in group),
            }
        )
    dump(folder / "summary.csv", result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    audit(p.parse_args().folder)
