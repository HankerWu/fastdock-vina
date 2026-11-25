#!/usr/bin/env python3
"""Check native output parity and mixed proposal accounting on 1IEP."""

import argparse
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory


FIXTURE = Path(__file__).resolve().parents[1] / "examples/basic"


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def search(binary, output, options):
    command = [
        str(binary), "--receptor", str(FIXTURE / "1iep_receptor.pdbqt"),
        "--ligand", str(FIXTURE / "1iep_ligand.pdbqt"),
        "--config", str(FIXTURE / "1iep_receptor.box.txt"), "--out", str(output),
        "--cpu", "1", "--seed", "42", "--exhaustiveness", "3", "--max_evals", "1800",
        "--num_modes", "3", "--energy_range", "100", "--verbosity", "0",
        *options,
    ]
    return subprocess.run(command, text=True, capture_output=True, timeout=90)


def verify_counters(log, mode, late):
    records = [
        dict(token.split("=", 1) for token in line.split()[2:])
        for line in log.splitlines() if line.startswith(("Search trajectory:", "Search counters:"))
    ]
    require(len(records) == 4, f"Expected three trajectory records and a summary: {log}")
    trajectories, summary = records[:3], records[3]
    require(summary.pop("rigid_reorientation") == mode, f"Wrong recorded mode: {log}")

    rotations = 0
    for index, record in enumerate(trajectories, start=1):
        counts = {key: int(value) for key, value in record.items()}
        require(counts["index"] == index, f"Out-of-order trajectory: {counts}")
        require(0 < counts["steps"] == counts["hunt_calls"] <= counts["hunt_evals"],
                f"Inconsistent hunt accounting: {counts}")
        require(0 < counts["second_calls"] <= counts["hunt_calls"]
                and counts["second_calls"] <= counts["second_evals"], f"Inconsistent second minimization: {counts}")
        require(counts["hunt_evals"] + counts["second_evals"] > 1800,
                f"Trajectory terminated before the requested score evaluations: {counts}")
        require(0 <= counts["rigid_reorientations"] <= counts["steps"],
                f"More global rotation proposals than steps: {counts}")
        rotations += counts["rigid_reorientations"]

    totals = {key: int(value) for key, value in summary.items()}
    require((rotations > 0) == (mode == "mixed"), f"Inconsistent mode and rotation count: {log}")
    require(0 < totals["explicit_inputs"] <= totals["explicit_calls"] <= totals["explicit_evals"],
            f"Inconsistent explicit refinement accounting: {totals}")
    require((totals["late_prescore"] > 0) == late, f"Unexpected Late prescore accounting: {totals}")
    print(f"{'Late' if late else 'Default'} {mode}: {rotations} global rotations, "
          f"{totals['explicit_calls']} explicit calls, {totals['explicit_evals']} explicit evaluations")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    baseline = args.baseline.resolve(strict=True)

    with TemporaryDirectory(prefix="vina-rigid-reorientation-") as work:
        directory = Path(work)

        def run(engine, name, flags):
            output = directory / f"{name}.pdbqt"
            result = search(engine, output, flags)
            require(result.returncode == 0, result.stderr or result.stdout)
            return output.read_bytes(), result.stdout

        for late in (False, True):
            options = ("--pose_selection", "late") if late else ()
            old = run(baseline, "old", options)
            default = run(binary, "default", options)
            native = run(binary, "native", (*options, "--rigid_reorientation", "native", "--search_counters"))
            mixed = run(binary, "mixed", (*options, "--rigid_reorientation", "mixed", "--search_counters"))
            replay = run(binary, "replay", (*options, "--rigid_reorientation", "mixed", "--search_counters"))
            require(old == default and native[0] == default[0], "Native output changed")
            require(mixed == replay and mixed[0] != native[0], "Mixed output must differ and be deterministic")
            verify_counters(native[1], "native", late)
            verify_counters(mixed[1], "mixed", late)

        rejected = search(binary, directory / "illegal.pdbqt", ("--local_only", "--rigid_reorientation", "mixed"))
        require(rejected.returncode != 0 and "mixed rigid_reorientation requires" in rejected.stderr,
                f"Unsupported local-only use did not fail clearly: {rejected.stderr}")


if __name__ == "__main__":
    main()
