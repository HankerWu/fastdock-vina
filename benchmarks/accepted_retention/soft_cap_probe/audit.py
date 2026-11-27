"""Validate archived diagnostic values against raw output and count pair inversions."""

import argparse
import csv
import hashlib
import statistics
import tarfile
from pathlib import Path


def read(path):
    with path.open() as file:
        return list(csv.DictReader(file))


def audit(folder):
    endpoints = read(folder / "endpoints.csv")
    cases = read(folder / "cases.csv")
    if len(cases) != 3 or len(endpoints) != 409:
        raise ValueError("Incomplete diagnostic panel")
    order = []
    with tarfile.open(folder / "runs.tar.xz") as tar:
        for case in cases:
            key = case["key"]
            values = [r for r in endpoints if r["key"] == key]
            by_step = {int(r["step"]): r for r in values}
            if len(by_step) != int(case["authentic_extra_joint_calls"]):
                raise ValueError("Incorrect evaluation count")
            for name, column in (
                ("parent", "parent_sha256"),
                ("candidate", "gated_sha256"),
                ("observer", "observer_sha256"),
            ):
                blob = tar.extractfile(f"{key}.{name}.pdbqt").read()
                if hashlib.sha256(blob).hexdigest() != case[column]:
                    raise ValueError("PDBQT hash differs")

            def counts(raw):
                prefixes = (
                    b"Native merge:",
                    b"Hunt chain:",
                    b"Search trajectory:",
                    b"Search counters:",
                    b"Deferred hunt:",
                )
                return [
                    l.lstrip(b"*")
                    for l in raw.splitlines()
                    if l.lstrip(b"*").startswith(prefixes)
                ]

            instrumented = counts(tar.extractfile(f"{key}.observer.stdout").read())
            native = counts(tar.extractfile(f"{key}.candidate.stdout").read())
            if (
                case["gated_sha256"] != case["observer_sha256"]
                or instrumented != native
                or len(native) != 5
            ):
                raise ValueError("Instrumented pool or counts differed")
            seen = set()
            for line in tar.extractfile(f"{key}.observer.stderr"):
                if not line.startswith(b"FDCAPS "):
                    continue
                _, index, soft, auth, dx2 = line.split()
                step = int(index)
                r = by_step[step]
                if (
                    step in seen
                    or float(dx2) != 0
                    or [float(soft), float(auth)]
                    != [float(r["soft"]), float(r["authentic"])]
                ):
                    raise ValueError("Raw evaluation differs")
                seen.add(step)
            if seen != by_step.keys():
                raise ValueError("Incomplete raw evaluations")
            pool = [
                r
                for r in values
                if min(float(r["nearest_merged_A"]), float(r["nearest_saved_A"])) >= 2
            ]
            if len(pool) != int(case["eligible_both"]):
                raise ValueError("Novel pool mismatch")
            ordered = sorted(pool, key=lambda r: (float(r["soft"]), int(r["step"])))
            hard_sorted = sorted(
                pool, key=lambda r: (float(r["authentic"]), int(r["step"]))
            )
            for metric, ranking in (("soft", ordered), ("authentic", hard_sorted)):
                oracle = [
                    i
                    for i, r in enumerate(ranking, 1)
                    if r["geometry_oracle"] == "True"
                ]
                if ranking[0]["step"] != case[f"{metric}_winner_step"] or oracle != [
                    int(case[f"oracle_{metric}_rank"])
                ]:
                    raise ValueError("Winner/rank differs")
            tied, discordant = 0, 0
            for i, x in enumerate(ordered):
                for y in ordered[:i]:
                    diff = (float(x["authentic"]) - float(y["authentic"])) * (
                        float(x["soft"]) - float(y["soft"])
                    )
                    tied += diff == 0
                    discordant += diff < 0
            pairs = len(pool) * (len(pool) - 1) // 2
            order.append(
                {
                    "key": key,
                    "eligible_states": len(pool),
                    "distinct_endpoint_pairs": pairs,
                    "tied": tied,
                    "inversions": discordant,
                    "median_delta_all": statistics.median(
                        float(x["delta"]) for x in values
                    ),
                    "max_delta_all": max(float(x["delta"]) for x in values),
                }
            )
    with (folder / "ordering.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=order[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(order)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results", type=Path)
    audit(p.parse_args().results)
