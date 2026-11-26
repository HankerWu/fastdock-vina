"""Check fixed-pool output equivalence and time matched command-line invocations."""

import argparse
import csv
import hashlib
import io
import tempfile
import resource
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from cli_checks import atom
from run import jobs, sample


PANELS = ("hcaii", "brd4", "t4l", "alkyl", "fxa", "fxa_pair", "bovine_ca", "trypsin")
VALUES = (
    "model",
    "inter",
    "intra",
    "slope",
    "window_delta",
    "half_delta",
    "uniform_delta",
    "effective_samples",
)
IDENTITY = ("receptor", "ligand", "seed", "arm", "slot", "model")


def verify(binary, base, frozen, samples, workers):
    for panel in PANELS:
        with (frozen / (panel + ".csv")).open() as stream:
            rows = list(csv.DictReader(stream))
        expected = {
            tuple(row[k] for k in IDENTITY): tuple(row[k] for k in VALUES)
            for row in rows
        }
        if len(expected) != len(rows):
            raise ValueError(f"{panel}: duplicate archived model")
        with ThreadPoolExecutor(workers) as pool:
            for output in pool.map(
                lambda job: sample(job, binary, samples), jobs(panel, base)
            ):
                for row in output:
                    key = tuple(row[k] for k in IDENTITY)
                    if expected.pop(key) != tuple(row[k] for k in VALUES):
                        raise ValueError(f"{panel}: changed printed quantity for {key}")
        if expected:
            raise ValueError(f"{panel}: missed {len(expected)} archived models")
        print(f"{panel}: {len(rows)} matched complete records", flush=True)


def timed_run(binary, job, samples, flags):
    _, receptor, ligand, poses, bonds, _ = job
    options = [
        binary,
        "--score_only",
        "--receptor",
        receptor,
        "--ligand",
        ligand,
        "--rigid_window_pool",
        poses,
        "--rigid_window_samples",
        str(samples),
        "--cpu",
        "1",
        "--verbosity",
        "0",
        *flags,
    ]
    if bonds:
        options.extend(["--ligand_bonds", bonds])
    start = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    stdout = subprocess.check_output(options)
    elapsed = time.perf_counter() - started
    finish = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = finish.ru_utime - start.ru_utime + finish.ru_stime - start.ru_stime
    return stdout[stdout.index(b"model,inter,intra,") :], cpu, elapsed


def verify_cutoff_edges(baseline, optimized):
    # Include the exact cutoff as well as its two nearest printed neighbours.
    with tempfile.TemporaryDirectory() as directory:
        rec, lig, pool = (
            Path(directory) / name for name in ("rec.pdbqt", "lig.pdbqt", "pool.pdbqt")
        )
        rec.write_text(atom(0, 0, 0))

        def ligand(distance):
            return "ROOT\n" + atom(distance, 0, 0) + "ENDROOT\nTORSDOF 0\n"

        for cutoff in (8.0, 7.73, 4.123):
            lig.write_text(ligand(cutoff))
            pool.write_text(
                "".join(
                    f"MODEL {i}\n" + ligand(cutoff + offset) + "ENDMDL\n"
                    for i, offset in enumerate((-0.001, 0, 0.001), 1)
                )
            )
            job = ((), rec, lig, pool, None, ())
            flags = ["--interaction_cutoff", str(cutoff)]
            before, _, _ = timed_run(baseline, job, 16, flags)
            after, _, _ = timed_run(optimized, job, 16, flags)
            inter = [
                float(r["inter"]) for r in csv.DictReader(io.StringIO(after.decode()))
            ]
            if before != after or inter[0] == 0 or inter[1:] != [0.0, 0.0]:
                raise ValueError(f"Incorrect cutoff boundary at {cutoff} A: {inter}")


def compare(base, baseline, optimized, samples, repeats, path):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            (
                "panel",
                "pool",
                "samples_per_scale",
                "mode",
                "repeat",
                "implementation",
                "models",
                "cpu_s",
                "wall_s",
                "csv_sha256",
            )
        )
        for panel in ("hcaii", "brd4", "trypsin", "fxa", "bovine_ca"):
            job = next(jobs(panel, base))
            for mode, flags in (
                ("native", []),
                ("custom_cutoff", ["--interaction_cutoff", "7.73"]),
            ):
                # At 7.73 A the cutoff straddles the last lookup bin.
                count = samples if mode == "native" else 128
                expected = None
                for rep in range(repeats if mode == "native" else 1):
                    versions = [("baseline", baseline), ("optimized", optimized)]
                    if rep % 2:
                        versions.reverse()
                    for name, binary in versions:
                        scores, cpu, wall = timed_run(binary, job, count, flags)
                        if expected is not None and scores != expected:
                            raise ValueError(
                                f"{job[3]} {mode}: complete model CSV differs"
                            )
                        expected = scores
                        writer.writerow(
                            (
                                panel,
                                job[3].name,
                                count,
                                mode,
                                rep,
                                name,
                                len(scores.splitlines()) - 1,
                                cpu,
                                wall,
                                hashlib.sha256(scores).hexdigest(),
                            )
                        )
                        stream.flush()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("base", type=Path)
    p.add_argument("optimized", type=Path)
    p.add_argument("--samples", type=int, default=2048)
    p.add_argument(
        "--frozen",
        type=Path,
        help="directory with previous panel CSVs at this sample count",
    )
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--baseline", type=Path)
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--timings", type=Path)
    args = p.parse_args()
    if bool(args.timings) != bool(args.baseline):
        p.error("provide both --baseline and --timings")
    if not args.frozen and not args.timings:
        p.error("specify --frozen or --timings")
    if min(args.repeats, args.jobs) < 1:
        p.error("repeats and jobs must be positive")
    if args.frozen:
        verify(
            args.optimized.resolve(), args.base, args.frozen, args.samples, args.jobs
        )
    if args.timings:
        verify_cutoff_edges(args.baseline.resolve(), args.optimized.resolve())
        compare(
            args.base,
            args.baseline.resolve(),
            args.optimized.resolve(),
            args.samples,
            args.repeats,
            args.timings,
        )
