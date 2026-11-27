"""Rescore every skipped accepted endpoint on three recovered native-pool failures."""

import argparse
import csv
import hashlib
import json
import os
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
COVERAGE = HERE.parent / "coverage"
sys.path.insert(0, str(COVERAGE))
from probe import counters, run_binary, table


def load(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(a):
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0)) if a.cpu is None else a.cpu})
    a.work.mkdir(parents=True, exist_ok=True)
    a.results.mkdir(parents=True, exist_ok=True)
    with tarfile.open(
        HERE.parents[1] / "cubic_grid/human_panel/prepared.tar.xz"
    ) as archive:
        archive.extractall(a.work, filter="data")
    meta = [r for r in load(COVERAGE / "results/opportunities.csv") if r["oracle_step"]]
    choices = load(COVERAGE / "results/choices.csv")
    frozen = {(r["key"], r["arm"]): r for r in load(COVERAGE / "results/blind.csv")}
    rows, cases = [], []
    raw = a.work / "raw"
    raw.mkdir(exist_ok=True)
    for spec in meta:
        key = spec["key"]
        known = {int(r["step"]): r for r in choices if r["key"] == key}
        logs, times, digests = {}, {}, {}
        for name in ("parent", "candidate", "observer"):
            dest = raw / f"{key}.{name}.pdbqt"
            p, elapsed = run_binary(a.work, spec, getattr(a, name), dest)
            logs[name], times[name], digests[name] = p, elapsed, sha(dest)
            (raw / f"{key}.{name}.stdout").write_bytes(p.stdout)
            (raw / f"{key}.{name}.stderr").write_bytes(p.stderr)
        if (
            digests["candidate"] != digests["observer"]
            or digests["candidate"] != frozen[key, "gated"]["sha256"]
            or digests["parent"] != frozen[key, "parent"]["sha256"]
        ):
            raise ValueError("Replayed pool hash differs")
        counted = counters(logs["candidate"].stdout)
        if len(counted) != 4 or counters(logs["observer"].stdout) != counted:
            raise ValueError("Observer changed production work")
        if [
            x for x in logs["candidate"].stdout.splitlines() if b"Deferred hunt:" in x
        ] != [
            x for x in logs["observer"].stdout.splitlines() if b"Deferred hunt:" in x
        ]:
            raise ValueError("Retention choice changed")
        parsed = {}
        for line in logs["observer"].stderr.splitlines():
            if not line.startswith(b"FDCAPS "):
                continue
            _, idx, e0, e1, dx2 = line.split()
            step = int(idx)
            r = known[step]
            soft, hard = float(e0), float(e1)
            if float(dx2) != 0 or soft != float(r["energy"]) or step in parsed:
                raise ValueError(
                    "Coordinate/soft-energy identity or uniqueness check failed"
                )
            obj = {
                "key": key,
                "step": step,
                "soft": soft,
                "authentic": hard,
                "delta": hard - soft,
                "rmsd_A": float(r["rmsd"]),
                "nearest_saved_A": float(r["nearest_saved"]),
                "nearest_merged_A": float(r["nearest_merged"]),
                "original_final_record": r["is_final_record"] == "True",
                "geometry_oracle": r["is_oracle"] == "True",
            }
            if hard < soft - 1e-12:
                raise ValueError(
                    "Larger positive-energy cap reduced this identical endpoint's energy"
                )
            parsed[step] = obj
        if parsed.keys() != known.keys():
            raise ValueError("Missing or duplicated skipped candidates")
        seq = list(parsed.values())
        novel = [r for r in seq if r["nearest_saved_A"] >= 2]
        # Whole archived pool permits checking whether a loser of the one-record competition was eligible at exit.
        eligible = [r for r in novel if r["nearest_merged_A"] >= 2]
        # PDBQT quantization changes each fixed-index RMSD by at most 0.0017321 A for two rounded coordinate arrays.
        if any(abs(r["nearest_merged_A"] - 2) <= 0.001733 for r in seq):
            raise ValueError(
                "Rounded final distance too close to classify; rerun with unrounded coordinates"
            )
        orig = next(r for r in seq if r["original_final_record"])
        target = next(r for r in seq if r["geometry_oracle"])
        summary = {
            "key": key,
            "skipped": len(seq),
            "novel_on_arrival": len(novel),
            "eligible_both": len(eligible),
            "authentic_extra_joint_calls": len(seq),
            "native_joint": int(frozen[key, "gated"]["joint"]),
            "max_movable_heavy_atom_displacement_A": 0,
            "parent_sha256": digests["parent"],
            "gated_sha256": digests["candidate"],
            "observer_sha256": digests["observer"],
            "parent_s": times["parent"],
            "gated_s": times["candidate"],
            "observer_s": times["observer"],
            "oracle_step": target["step"],
            "oracle_rmsd_A": target["rmsd_A"],
            "oracle_soft": target["soft"],
            "oracle_authentic": target["authentic"],
        }
        for metric in ("soft", "authentic"):
            record = min(novel, key=lambda r: (r[metric], r["step"]))
            winner = min(eligible, key=lambda r: (r[metric], r["step"]))
            if metric == "soft" and record != orig:
                raise ValueError("Failed to reproduce online selected reserve")
            summary.update(
                {
                    f"{metric}_record_step": record["step"],
                    f"{metric}_record_survives": record in eligible,
                    f"{metric}_winner_step": winner["step"],
                    f"{metric}_winner_rmsd_A": winner["rmsd_A"],
                    f"{metric}_winner_energy": winner[metric],
                    f"oracle_{metric}_rank": sorted(
                        eligible, key=lambda r: (r[metric], r["step"])
                    ).index(target)
                    + 1,
                    f"oracle_{metric}_gap": target[metric] - winner[metric],
                }
            )
        rows.extend(seq)
        cases.append(summary)
    table(a.results / "endpoints.csv", rows)
    table(a.results / "cases.csv", cases)
    (a.results / "binary_sha256.json").write_text(
        json.dumps(
            {n: sha(getattr(a, n)) for n in ("parent", "candidate", "observer")},
            indent=2,
        )
        + "\n"
    )
    with tarfile.open(a.results / "runs.tar.xz", "w:xz") as archive:
        for f in sorted(raw.iterdir()):
            archive.add(f, arcname=f.name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "results", "parent", "candidate", "observer"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--cpu", type=int)
    run(parser.parse_args())
