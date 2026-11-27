"""Trace native single-trajectory runs on all 66 fixed six-family inputs."""

import argparse
import csv
import hashlib
import io
import json
import os
import subprocess
import sys
import tarfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
sys.path.insert(0, str(BENCH / "paired_crossdocking"))
from run import alignment

PANEL = BENCH / "cubic_grid/human_panel"
STAGES = (
    "start",
    "proposal",
    "hunt_accepted",
    "hunt_rejected",
    "saved_retained",
    "saved_dropped",
    "mc_exit",
    "merged",
    "explicit_input",
    "explicit_output",
)


def table(path, rows):
    with path.open("w") as f:
        writer = csv.DictWriter(f, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def command(work, row, binary, dest):
    q, r = row["query"], row["receiver"]
    rm = json.loads((work / r / "preparation.json").read_text())
    args = [
        binary,
        "--ligand",
        work / q / "ligand.pdbqt",
        "--receptor",
        work / r / "receptor.pdbqt",
        "--seed",
        str(row["seed"]),
        "--max_evals",
        "12000",
        "--num_modes",
        "9",
        "--cpu",
        "1",
        "--exhaustiveness",
        "1",
        "--energy_range",
        "100",
        "--search_counters",
        "--out",
        dest,
    ]
    for axis, x in zip("xyz", rm["center"]):
        args += ["--center_" + axis, str(x), "--size_" + axis, "24"]
    return args


def run_binary(work, row, binary, dest):
    start = time.perf_counter()
    p = subprocess.run(
        command(work, row, binary, dest), capture_output=True, check=True
    )
    return p, time.perf_counter() - start


def references(work, row):
    q, r = [
        json.loads((work / row[k] / "preparation.json").read_text())
        for k in ("query", "receiver")
    ]
    rot, shift, _, _ = alignment(q, r)
    return q["serial"], np.asarray(q["mapped"]) @ rot + shift


def rmsd(x, refs):
    return np.sqrt(np.square(refs - x).sum(axis=2).mean(axis=1)).min()


def mode_distances(data, serial, refs):
    result = []
    for block in data.decode().split("ENDMDL")[:-1]:
        pos = {
            int(s[6:11]): [float(s[a : a + 8]) for a in (30, 38, 46)]
            for s in block.splitlines()
            if s.startswith(("ATOM  ", "HETATM"))
        }
        result.append(rmsd(np.array([pos[i] for i in serial]), refs))
    return result


def counters(log):
    return [
        s
        for s in log.splitlines()
        if s.startswith(
            (
                b"Native merge:",
                b"Hunt chain:",
                b"Search trajectory:",
                b"Search counters:",
            )
        )
    ]


def measure(work, row, parent, tracer):
    key = "_".join(row[k] for k in ("target", "query", "receiver", "seed"))
    dest = work / "runs" / key
    plain, plain_s = run_binary(work, row, parent, dest.with_suffix(".pdbqt"))
    trace, trace_s = run_binary(work, row, tracer, dest.with_suffix(".traced"))
    data = dest.with_suffix(".pdbqt").read_bytes()
    if (
        data != dest.with_suffix(".traced").read_bytes()
        or len(counters(plain.stdout)) != 4
        or counters(plain.stdout) != counters(trace.stdout)
    ):
        raise ValueError("Passive parity failed: " + key)
    dest.with_suffix(".log").write_bytes(plain.stdout + plain.stderr)
    serial, refs = references(work, row)
    positions, stages, steps, values = [], [], [], []
    order = None
    for line in trace.stderr.splitlines():
        if not line.startswith(b"FDTRACE "):
            continue
        _, stage, idx, xyz = line.split(maxsplit=3)
        # The trace map is sorted by PDBQT serial. Each atom has id,x,y,z.
        atoms = np.fromstring(xyz.replace(b" ", b","), sep=",").reshape(-1, 4)
        if order is None:
            ids = atoms[:, 0].astype(int).tolist()
            if set(ids) != set(serial) or len(ids) != len(serial):
                raise ValueError("Atom map mismatch")
            order = [ids.index(i) for i in serial]
        if not np.array_equal(atoms[order, 0], serial):
            raise ValueError("Trace changed atom order")
        xyz = atoms[order, 1:]
        positions.append(np.rint(1000 * xyz).astype(np.int32))
        stages.append(STAGES.index(stage.decode()))
        steps.append(int(idx))
        values.append(rmsd(xyz, refs))
    # Store coordinates losslessly at PDBQT precision; retain full precision RMSD for threshold cases.
    np.savez(
        dest.with_suffix(".npz"),
        stage=np.array(stages, np.uint8),
        step=np.array(steps, np.int32),
        xyz_milli=np.stack(positions),
        rmsd=np.array(values, np.float64),
    )
    dest.with_suffix(".traced").unlink()
    output = mode_distances(data, serial, refs)
    case = {k: row[k] for k in ("target", "query", "receiver", "kind", "seed")}
    case.update(
        key=key,
        events=len(steps),
        modes=len(output),
        top1=output[0],
        top5=min(output[:5]),
        oracle=min(output),
        native_s=plain_s,
        trace_s=trace_s,
        native_sha256=hashlib.sha256(data).hexdigest(),
        trace_sha256=hashlib.sha256(trace.stderr).hexdigest(),
    )
    return case


def archive(work, out, group, cases):
    # Normalized metadata for reproducible compressed archives.
    with tarfile.open(out / (group + ".tar.xz"), "w:xz", preset=4) as t:
        for case in cases:
            for ext in (".npz", ".log", ".pdbqt"):
                p = work / "runs" / (case["key"] + ext)
                data = p.read_bytes()
                info = tarfile.TarInfo(p.name)
                info.size = len(data)
                t.addfile(info, io.BytesIO(data))


def main(args):
    os.sched_setaffinity(0, {args.cpu})
    for path in (args.work / "runs", args.output):
        path.mkdir(parents=True, exist_ok=True)
    with tarfile.open(PANEL / "prepared.tar.xz") as t:
        t.extractall(args.work, filter="data")
    rows = list(csv.DictReader((PANEL / "cells.csv").open()))
    assert len(rows) == 66
    cases = []
    for r in rows:
        cases.append(measure(args.work, r, args.parent, args.trace))
    table(args.output / "cells.csv", cases)
    for group in sorted({r["target"] for r in cases}):
        archive(
            args.work, args.output, group, [r for r in cases if r["target"] == group]
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "output", "parent", "trace"):
        p.add_argument(name, type=Path)
    p.add_argument("--cpu", type=int, required=True)
    main(p.parse_args())
