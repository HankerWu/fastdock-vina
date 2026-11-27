"""Verify coordinates and identify when the native single trajectory loses a 2 Å hit."""

import argparse
import csv
import hashlib
import re
import tarfile
from pathlib import Path

import numpy as np
from probe import PANEL, STAGES, counters, mode_distances, references, table


def measure(work, root):
    rows = list(csv.DictReader((root / "cells.csv").open()))
    if len(rows) != 66 or len({r["key"] for r in rows}) != 66:
        raise ValueError("Unexpected experimental cells")
    with tarfile.open(PANEL / "prepared.tar.xz") as t:
        t.extractall(work, filter="data")
    for group in sorted({r["target"] for r in rows}):
        with tarfile.open(root / (group + ".tar.xz")) as t:
            t.extractall(work / "runs", filter="data")
    results = []
    for row in rows:
        prefix = work / "runs" / row["key"]
        data = prefix.with_suffix(".pdbqt").read_bytes()
        archive = np.load(prefix.with_suffix(".npz"), allow_pickle=False)
        names = np.asarray(STAGES)[archive["stage"]]
        steps, errors = archive["step"], archive["rmsd"]
        serial, refs = references(work, row)
        xyz = archive["xyz_milli"].astype(float) / 1000
        checked = []
        for block in np.array_split(xyz, max(1, len(xyz) // 48)):
            # Rows are raw deposited heavy-atom order, not ligand-fit coordinates.
            dd = block[:, None, :, :] - refs
            checked.extend(np.sqrt(np.square(dd).sum(axis=3).mean(axis=2)).min(axis=1))
        if max(abs(checked - errors)) > np.sqrt(3) * 0.0005 + 1e-10:
            raise ValueError("Coordinate-ledger disagreement")
        if not np.array_equal(np.array(checked) <= 2, errors <= 2):
            raise ValueError("Coordinate precision changes classification")
        score = mode_distances(data, serial, refs)
        if (
            hashlib.sha256(data).hexdigest() != row["native_sha256"]
            or abs(min(score) - float(row["oracle"])) > 1e-10
        ):
            raise ValueError("Final pool differs from capture")
        counts = dict(
            re.findall(
                rb"(\w+)=(\d+)",
                b" ".join(counters(prefix.with_suffix(".log").read_bytes())),
            )
        )
        accepted = names == "hunt_accepted"
        submitted_steps = steps[np.isin(names, ["saved_retained", "saved_dropped"])]
        submitted = accepted & np.isin(steps, submitted_steps)
        skipped = accepted & ~submitted
        if (
            int(counts[b"accepted"]) != sum(accepted)
            or int(counts[b"second_calls"]) != len(submitted_steps)
            or int(counts[b"hunt_calls"]) != sum(names == "proposal")
        ):
            raise ValueError("Trace not exhaustive")
        entry = {
            k: row[k] for k in ("key", "target", "query", "receiver", "kind", "seed")
        }
        groups = {
            "proposal": names == "proposal",
            "hunt_rejected": names == "hunt_rejected",
            "accepted": accepted,
            "submitted": submitted,
            "skipped": skipped,
            "second": np.isin(names, ["saved_retained", "saved_dropped"]),
            "saved": names == "saved_retained",
            "exit": names == "mc_exit",
            "merged": names == "merged",
            "explicit": names == "explicit_output",
        }
        for label, mask in groups.items():
            entry[label] = float(min(errors[mask])) if mask.any() else ""
        entry.update(
            n_accepted=int(accepted.sum()),
            n_submitted=int(submitted.sum()),
            n_skipped=int(skipped.sum()),
            top1=score[0],
            top5=min(score[:5]),
            oracle=min(score),
        )
        eligible = skipped & (errors <= 2)
        idx = np.where(eligible)[0]
        entry["oracle_step"] = (
            int(steps[idx[np.argmin(errors[idx])]])
            if len(idx) and min(score) > 2
            else ""
        )
        results.append(entry)
    return results


def summary(rows):
    def hit(r, k):
        return r[k] != "" and float(r[k]) <= 2

    for group in ["all", *sorted({r["target"] for r in rows})]:
        used = rows if group == "all" else [r for r in rows if r["target"] == group]
        entry = {"target": group, "cells": len(used)}
        for field in (
            "proposal",
            "hunt_rejected",
            "accepted",
            "submitted",
            "skipped",
            "second",
            "saved",
            "exit",
            "merged",
            "explicit",
            "top1",
            "top5",
            "oracle",
        ):
            entry[field] = sum(hit(r, field) for r in used)
        entry.update(
            opportunity=sum(r["oracle_step"] != "" for r in used),
            skipped_only=sum(
                r["oracle_step"] != "" and not hit(r, "submitted") for r in used
            ),
            endpoint_absent=sum(
                not hit(r, "accepted")
                and not hit(r, "hunt_rejected")
                and not hit(r, "oracle")
                for r in used
            ),
        )
        yield entry


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("results", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rows = measure(args.work, args.results)
    table(args.output / "opportunities.csv", rows)
    table(args.output / "coverage.csv", list(summary(rows)))
