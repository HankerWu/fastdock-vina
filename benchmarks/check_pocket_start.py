#!/usr/bin/env python3
"""Check optional initial pocket placement on real 1IEP inputs."""

import argparse
from math import dist
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory


FIXTURE = Path(__file__).resolve().parents[1] / "examples/basic"


def require(condition, reason):
    if not condition:
        raise RuntimeError(reason)


def telemetry(log, prefix):
    return {
        int(values.pop("index")): values
        for line in log.splitlines() if line.startswith(prefix)
        if (values := dict(token.split("=", 1) for token in line[len(prefix):].split()))
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path)
    args = parser.parse_args()
    binary, baseline = (item.resolve(strict=True) for item in (args.binary, args.baseline))
    box = FIXTURE / "1iep_receptor.box.txt"
    box_values = {key.strip(): float(value) for key, value in (
        line.split("=") for line in box.read_text().splitlines() if line.strip()
    )}
    center = tuple(box_values[f"center_{axis}"] for axis in "xyz")
    with TemporaryDirectory(prefix="vina-pocket-start-") as work:
        scratch = Path(work)

        def search(engine, name, flags, cpu=1, insist=True):
            output = scratch / f"{name}.pdbqt"
            command = [str(engine), "--receptor", str(FIXTURE / "1iep_receptor.pdbqt"),
                       "--ligand", str(FIXTURE / "1iep_ligand.pdbqt"), "--config", str(box),
                       "--out", str(output), "--cpu", str(cpu), "--seed", "42", "--max_evals", "3000",
                       "--num_modes", "3", "--exhaustiveness", "4", "--verbosity", "0", *flags]
            result = subprocess.run(command, text=True, capture_output=True, timeout=90)
            if insist:
                require(result.returncode == 0, result.stderr or result.stdout)
                return output.read_bytes(), result.stdout
            return result

        old = search(baseline, "old", [])
        default = search(binary, "default", [])
        require(old == default, "The default 1IEP PDBQT or stdout changed")
        logged = search(binary, "logged", ["--search_counters"])
        require(logged[0] == default[0], "Observing pocket starters changed a native pool")
        flags = ["--search_counters", "--pocket_start_radius", "3"]
        mixed = search(binary, "pocket", flags)
        threaded = search(binary, "threaded", flags, cpu=2)
        require(mixed == threaded, "Pocket starts depend on worker scheduling")
        native_starts = telemetry(logged[1], "Pocket start:")
        mixed_starts = telemetry(mixed[1], "Pocket start:")
        native_steps = telemetry(logged[1], "Search trajectory:")
        mixed_steps = telemetry(mixed[1], "Search trajectory:")
        require(list(mixed_starts) == [1, 2, 3, 4], f"Missing trajectories: {mixed_starts}")
        require(all(row["anchored"] == "0" for row in native_starts.values()), "Native trajectory anchored")
        for index, row in mixed_starts.items():
            anchored = index % 2 == 0
            require(int(row["anchored"]) == anchored, "Unexpected anchoring schedule")
            if anchored:
                centroid = tuple(float(row[f"c{axis}"]) for axis in "xyz")
                require(dist(centroid, center) <= 3 + 1e-10, "Anchored centroid lies outside the sphere")
            else:
                require(row == native_starts[index] and mixed_steps[index] == native_steps[index],
                        "An untouched trajectory changed its starting coordinates or search work")

        for invalid in ("-1", "0", "nan", "inf", "11"):
            error = search(binary, f"invalid-{invalid}", ["--pocket_start_radius", invalid], insist=False)
            require(error.returncode != 0 and "pocket_start_radius" in error.stderr,
                    f"Invalid radius accepted: {invalid}")
        wrong_mode = search(binary, "local", ["--pocket_start_radius", "3", "--local_only"], insist=False)
        require(wrong_mode.returncode != 0 and "pocket_start_radius" in wrong_mode.stderr,
                "Local-only docking accepted a pocket starting radius")
    print("1IEP native byte parity, pocket geometry, untouched trajectory counters and CLI validation passed.")


if __name__ == "__main__":
    main()
