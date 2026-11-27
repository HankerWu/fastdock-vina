"""Replay transfer archive updates, fixed-slot ablations, or a common 11,400 cap."""

import argparse
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import run_panel

sys.path.insert(0, str(HERE.parents[1] / "development"))
import replay


def arguments(source, meta, q, r, seed, cap):
    common = [
        "--receptor",
        str(source / r / "receptor.pdbqt"),
        "--ligand",
        str(source / q / "ligand.pdbqt"),
        "--cpu",
        "1",
        "--exhaustiveness",
        "4",
        "--max_evals",
        str(cap),
        "--num_modes",
        "9",
        "--min_rmsd",
        "1.0",
        "--energy_range",
        "100",
        "--seed",
        str(seed),
    ]
    for axis, value in zip("xyz", meta[r]["center"]):
        common += ["--center_" + axis, str(value), "--size_" + axis, "22.5"]
    return common


def refinement_slots(log, output, key, query, receiver):
    sys.path.insert(0, str(HERE.parents[2] / "pair_gradient/holdout"))
    from run_holdout import TRACE, stage_data

    def geometry(block):
        atoms = run_panel.paired.atoms(block.splitlines())
        return tuple(tuple(atoms[i][0]) for i in query["serial"])

    ranks = {
        geometry(block): i for i, block in enumerate(output.split("ENDMDL")[:-1], 1)
    }
    mapped = {
        int(index): ranks.get(geometry(block.split("PDBQT\n")[1]), "")
        for stage, index, block in TRACE.findall(log)
        if stage == "after"
    }
    if sorted(r for r in mapped.values() if r != "") != list(range(1, len(ranks) + 1)):
        raise ValueError("Ambiguous slot mapping: " + key)
    rows = stage_data(log, key, query, receiver)
    for row in rows:
        row["final_rank"] = (
            mapped[row["input_index"]] if row["stage"] == "after" else ""
        )
    return rows


def run(args):
    source, dest = args.panel.resolve(), args.output.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    raw = dest / "runs"
    raw.mkdir(exist_ok=True)
    os.sched_setaffinity(0, {args.cpu})
    keys = [
        (q, r, s)
        for q, r in [
            ("3P6D", "3P6D"),
            ("3P6E", "3P6D"),
            ("3P6D", "3P6E"),
            ("3P6E", "3P6E"),
            ("3H82", "3H82"),
            ("3F1O", "3H82"),
            ("3H82", "3F1O"),
            ("3F1O", "3F1O"),
        ]
        for s in run_panel.SEEDS
    ]
    meta = {
        q: json.loads((source / q / "preparation.json").read_text()) for q, _, _ in keys
    }
    for q in meta:
        if not (dest / q).exists():
            (dest / q).symlink_to(source / q, target_is_directory=True)
    arms = ("parent", "candidate") if args.mode != "factorial" else ("save", "merge")
    bins = {k: getattr(args, k).resolve() for k in arms}
    (dest / "binaries.json").write_text(
        json.dumps({k: run_panel.sha(v) for k, v in bins.items()}, indent=2) + "\n"
    )
    rows, modes, streams, changes = [], [], [], []
    cap = 11400 if args.mode == "reduced" else 12000
    for i, (q, r, seed) in enumerate(keys):
        common = arguments(source, meta, q, r, seed, cap)
        for arm in arms[:: (-1 if i % 2 else 1)]:
            key = f"{q}_{r}_{seed}_{arm}"
            output, log = raw / (key + ".pdbqt"), raw / (key + ".log")
            row = dict(
                key=key,
                query=q,
                receiver=r,
                seed=seed,
                cap=cap,
                arm=arm,
                target=meta[q]["target"],
            )
            if args.mode == "refinement":
                with log.open("w") as handle:
                    subprocess.run(
                        [bins[arm], *common, "--out", str(output), "--search_counters"],
                        stdout=handle,
                        stderr=subprocess.STDOUT,
                        check=True,
                    )
                if output.read_bytes() != (source / "runs" / output.name).read_bytes():
                    raise ValueError(key + " pose mismatch")
                traced = run_panel.sum_counts(log.read_text())
                if traced != run_panel.sum_counts(
                    (source / "runs" / (key + ".counts")).read_text()
                ):
                    raise ValueError(key + " count mismatch")
                streams.extend(
                    refinement_slots(
                        log.read_text(), output.read_text(), key, meta[q], meta[r]
                    )
                )
                row.update(traced)
            elif args.mode == "trace":
                trace = raw / (key + ".trace")
                with log.open("w") as handle:
                    subprocess.run(
                        [bins[arm], *common, "--out", str(output), "--search_counters"],
                        stdout=handle,
                        stderr=subprocess.STDOUT,
                        check=True,
                        env={**os.environ, "FD_ARCHIVE_TRACE": str(trace)},
                    )
                if output.read_bytes() != (source / "runs" / output.name).read_bytes():
                    raise ValueError(key + " pose mismatch")
                traced = run_panel.sum_counts(log.read_text())
                if traced != run_panel.sum_counts(
                    (source / "runs" / (key + ".counts")).read_text()
                ):
                    raise ValueError(key + " count mismatch")
                packed = Path(str(trace) + ".gz")
                with trace.open("rb") as i_stream, gzip.open(packed, "wb") as o_stream:
                    shutil.copyfileobj(i_stream, o_stream)
                trace.unlink()
                row.update(traced)
                st, ch = replay.scan(
                    packed,
                    meta[q],
                    meta[r],
                    key,
                    "parent" if arm == "parent" else "fixed",
                )
                streams.extend(st)
                changes.extend(ch)
            else:
                wall, cpu, rss = run_panel.cli(bins[arm], common, output, log)
                counted = raw / (key + ".replay")
                cntlog = raw / (key + ".counts")
                run_panel.cli(bins[arm], common, counted, cntlog, counters=True)
                if output.read_bytes() != counted.read_bytes():
                    raise ValueError(key + " counter replay mismatch")
                counted.unlink()
                row.update(
                    wall_s=wall,
                    cpu_s=cpu,
                    maxrss_KiB=rss,
                    **run_panel.sum_counts(cntlog.read_text()),
                )
            modes.extend(run_panel.paired.inspect(raw, meta[q], meta[r], row))
            row["sha256"] = run_panel.sha(output)
            rows.append(row)
            print(key, flush=True)
    run_panel.paired.write(dest / "cells.csv", rows)
    run_panel.paired.write(dest / "modes.csv", modes)
    if streams:
        if args.mode == "refinement":
            streams.sort(key=lambda row: row["key"])
        run_panel.paired.write(
            dest / ("slots.csv" if args.mode == "refinement" else "streams.csv"),
            streams,
        )
    if changes:
        run_panel.paired.write(dest / "interventions.csv", changes)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("panel", "output"):
        p.add_argument(name, type=Path)
    p.add_argument("mode", choices=("trace", "reduced", "factorial", "refinement"))
    p.add_argument("--cpu", type=int, required=True)
    for arm in ("parent", "candidate", "save", "merge"):
        p.add_argument("--" + arm, type=Path)
    run(p.parse_args())
