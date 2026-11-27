"""Replay both bounded archive policies on identical full-precision candidate streams."""

import argparse
import csv
import gzip
import itertools
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "rigid_flip"))
import holdout


def distance(x, y):
    return float(np.sqrt(np.square(x - y).sum(axis=1).mean()))


def update(pool, p, radius, cap, xyz, policy):
    ds = [distance(xyz[a], xyz[p]) for a in pool]
    nearby = [a for a, d in zip(pool, ds) if d < radius]
    if policy == "nearest" and nearby:
        i = int(np.argmin(ds))
        if p[0] >= pool[i][0]:
            return pool.copy()
        selected = pool.copy()
        selected[i] = p
    else:
        if any(a[0] <= p[0] for a in nearby):
            return pool.copy()
        selected = [a for a in pool if a not in nearby] + [p]
    return sorted(selected, key=lambda k: k[0])[:cap]


def signatures(tokens):
    return [tuple(map(float, t.split(","))) for t in tokens.split()]


def scan(trace, q, r, key, arm):
    rotation, shift, _, _ = holdout.paired.base.alignment(q, r)
    reference = np.array(q["mapped"]) @ rotation + shift
    xyz, rms = {}, {}
    streams, changed = [], []
    pool = alt = []
    current_stage = None

    def finish():
        if not streams:
            return
        row = streams[-1]
        row.update(
            end_size=len(pool),
            counterfactual_size=len(alt),
            end_oracle_A=min(rms[a] for a in pool),
            counterfactual_oracle_A=min(rms[a] for a in alt),
            same_end=int(pool == alt),
        )

    opener = gzip.open if trace.exists() else open
    trace = trace if trace.exists() else trace.with_suffix("")
    with opener(trace, "rt") as src:
        for line in src:
            if line.startswith("atoms "):
                serial = list(map(int, line.split()[1:]))
                if sorted(serial) != sorted(q["serial"]):
                    raise ValueError("Unexpected heavy-atom ordering")
                order = [serial.index(i) for i in q["serial"]]
                continue
            trial, before, after = line.rstrip("\n").split("\t")
            stage, radius, cap, energy, coordinates = trial.split(" ", 4)
            radius, cap, energy = float(radius), int(cap), float(energy)
            positions = np.fromstring(coordinates, sep=" ").reshape(-1, 3)
            identifier = (energy, *positions[0])
            if identifier in xyz and not np.array_equal(xyz[identifier], positions):
                raise ValueError("Pose key collision")
            xyz[identifier] = positions
            rms[identifier] = np.sqrt(
                np.square(reference - positions[order]).sum(axis=2).mean(axis=1)
            ).min()
            before, after = signatures(before), signatures(after)
            if not before:
                finish()
                pool, alt = [], []
                streams.append(
                    dict(
                        key=key,
                        arm=arm,
                        stream=len(streams) + 1,
                        stage=stage,
                        radius=radius,
                        arrivals=0,
                        multi_neighbor=0,
                        immediate_difference=0,
                        invalid_steps=0,
                        max_violations=0,
                        minimum_gap_A=np.inf,
                        near_native_inputs=0,
                    )
                )
                current_stage = stage
            if pool != before or current_stage != stage:
                raise ValueError(f"Failed to reproduce pre-state in {key}")
            row = streams[-1]
            row["arrivals"] += 1
            row["near_native_inputs"] += int(rms[identifier] <= 2)
            nearby = [
                (a, distance(xyz[a], positions))
                for a in before
                if distance(xyz[a], positions) < radius
            ]
            row["multi_neighbor"] += int(len(nearby) > 1)
            left = update(before, identifier, radius, cap, xyz, "nearest")
            right = update(before, identifier, radius, cap, xyz, "all")
            if left != right:
                row["immediate_difference"] += 1
                changed.append(
                    dict(
                        key=key,
                        arm=arm,
                        stage=stage,
                        stream=row["stream"],
                        arrival=row["arrivals"],
                        radius_A=radius,
                        capacity=cap,
                        old_occupancy=len(before),
                        neighbors=len(nearby),
                        newcomer_energy=energy,
                        newcomer_rmsd_A=rms[identifier],
                        neighbor_energy_distance_rmsd=json.dumps(
                            [(a[0], d, rms[a]) for a, d in nearby]
                        ),
                        nearest_size=len(left),
                        all_size=len(right),
                        nearest_best_A=min(rms[a] for a in left),
                        all_best_A=min(rms[a] for a in right),
                    )
                )
            pool = left if arm == "parent" else right
            if pool != after:
                raise ValueError(
                    f"Full precision post-state mismatch in {key}: {pool} {after}"
                )
            alt = update(
                alt,
                identifier,
                radius,
                cap,
                xyz,
                "all" if arm == "parent" else "nearest",
            )
            sep = [distance(xyz[a], xyz[b]) for a, b in itertools.combinations(pool, 2)]
            violations = sum(d < radius for d in sep)
            row["max_violations"] = max(row["max_violations"], violations)
            row["invalid_steps"] += int(bool(violations))
            if sep:
                row["minimum_gap_A"] = min(row["minimum_gap_A"], min(sep))
        finish()
    if len(streams) != 6 or [r["stage"] for r in streams] != ["save"] * 4 + [
        "merge",
        "prefine",
    ]:
        raise ValueError(
            "Expected four search streams, one merge, one pre-refinement stream"
        )
    return streams, changed


def main(args):
    all_streams, all_changes = [], []
    for row in csv.DictReader((args.results / "cells.csv").open()):
        load = lambda k: json.loads(
            (args.inputs / row[k] / "preparation.json").read_text()
        )
        s, c = scan(
            args.results / "raw" / (row["key"] + ".trace.gz"),
            load("query"),
            load("receiver"),
            row["key"],
            row["arm"],
        )
        all_streams.extend(s)
        all_changes.extend(c)
    holdout.paired.write(args.results / "streams.csv", all_streams)
    holdout.paired.write(args.results / "interventions.csv", all_changes)
    print(
        f"{len(all_streams)} streams, {sum(s['arrivals'] for s in all_streams)} verified updates"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", type=Path)
    p.add_argument("results", type=Path)
    main(p.parse_args())
