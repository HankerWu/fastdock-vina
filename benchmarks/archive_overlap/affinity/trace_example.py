"""Locate the first unequal archive update in two paired, full-precision C++ recordings."""

import argparse
import collections
import json
import tarfile
from pathlib import Path

import numpy as np


def summarize(path):
    with tarfile.open(path) as data:
        runs = [
            data.extractfile(arm + ".trace").read().decode().splitlines()
            for arm in ("parent", "candidate")
        ]
    positions = {}
    stage_counts = collections.Counter()
    trajectory = 0
    for left, right in zip(*runs):
        if left.startswith("atoms"):
            assert left == right
            trajectory += 1
            continue
        candidate, before, after = left.split("\t")
        stage, radius, capacity, energy, coords = candidate.split(" ", 4)
        coords = np.fromstring(coords, sep=" ").reshape(-1, 3)
        positions[(float(energy), *coords[0])] = coords
        stage_counts[trajectory, stage] += 1
        if left == right:
            continue
        new_candidate, new_before, new_after = right.split("\t")
        assert candidate == new_candidate and before == new_before
        neighbors = []
        for signature in before.split():
            key = tuple(map(float, signature.split(",")))
            distance = np.sqrt(np.square(positions[key] - coords).sum(axis=1).mean())
            if distance < float(radius):
                neighbors.append({"energy": key[0], "distance_A": float(distance)})
        return {
            "trajectory": trajectory,
            "stage": stage,
            "arrival": stage_counts[trajectory, stage],
            "radius_A": float(radius),
            "capacity": int(capacity),
            "energy": float(energy),
            "neighbors": neighbors,
            "before": len(before.split()),
            "after_parent": len(after.split()),
            "after_candidate": len(new_after.split()),
            "total_arrivals": [
                dict(
                    collections.Counter(
                        line.split()[0] for line in run if not line.startswith("atoms")
                    )
                )
                for run in runs
            ],
        }
    raise ValueError(f"No differing archive state in {path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archives", type=Path, nargs="+")
    print(
        json.dumps(
            {path.name: summarize(path) for path in p.parse_args().archives}, indent=2
        )
    )
