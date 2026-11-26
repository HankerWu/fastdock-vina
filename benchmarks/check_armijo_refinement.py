#!/usr/bin/env python3
"""Check Armijo opt-in behavior, default parity and composition on a real Vina case."""

import argparse
import math
from pathlib import Path
import re
import subprocess
import tempfile

from check_late_pose_selection import COUNTS, ROOT, dock, parse, require


LOCAL_REFINEMENT = re.compile(r"^Local refinement: optimizer=(\w+) gradient_evals=(\d+) objective_kcal_mol=(\S+)$", re.M)


def coordinates(text):
    return tuple(tuple(map(float, (line[30:38], line[38:46], line[46:54])))
                 for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM")))


def check_local_budget(binary, workdir):
    fixture = ROOT / "examples/basic"
    ligand = fixture / "1iep_ligand.pdbqt"

    def optimize(mode, cap, input_pose=ligand, expected_error=None):
        output = workdir / f"{input_pose.stem}-{mode}-{cap}.pdbqt"
        command = [str(binary), "--receptor", str(fixture / "1iep_receptor.pdbqt"),
                   "--ligand", str(input_pose), "--autobox", "--local_only", "--cpu", "1",
                   "--seed", "42", "--verbosity", "0", "--refinement_optimizer", mode,
                   "--max_evals", str(cap), "--out", str(output)]
        result = subprocess.run(command, text=True, capture_output=True, timeout=90)
        if expected_error:
            require(result.returncode != 0 and expected_error in result.stderr and not output.exists(),
                    f"Unexpected negative budget handling: {result.stderr}")
            return
        require(result.returncode == 0, result.stderr or result.stdout)
        (reported_mode, count, objective), = LOCAL_REFINEMENT.findall(result.stdout)
        require(reported_mode == mode, f"Expected {mode}, got {reported_mode}")
        require(math.isfinite(float(objective)), "Local objective must be finite")
        return int(count), float(objective), coordinates(output.read_text())

    optimize("geometry", -1, expected_error="non-negative max_evals")
    initial = coordinates(ligand.read_text())
    any_descent = False
    for mode in ("armijo", "geometry", "rigid_metric"):
        one = optimize(mode, 1)
        short = optimize(mode, 7)
        unlimited = optimize(mode, 0)
        require(one[0] == 1 and one[2] == initial, f"{mode} changed coordinates after one gradient")
        require(1 < short[0] <= 7 < unlimited[0], f"{mode} did not observe the finite gradient cap")
        require(unlimited[1] <= short[1] <= one[1], f"{mode} did not preserve the accepted objective")
        any_descent |= unlimited[1] < one[1]
    require(any_descent, "The real fixture did not exercise an accepted trial step")

    carbon = workdir / "one-carbon.pdbqt"
    atom = next(line for line in ligand.read_text().splitlines()
                if line.startswith("ATOM  ") and line.split()[-1] == "C")
    carbon.write_text(f"ROOT\n{atom}\nENDROOT\nTORSDOF 0\n")
    baseline = optimize("armijo", 7, carbon)
    for mode in ("geometry", "rigid_metric"):
        require(optimize(mode, 7, carbon) == baseline, f"{mode} changed Armijo on a single atom")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("pre_optimizer", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    baseline = args.pre_optimizer.resolve(strict=True)
    with tempfile.TemporaryDirectory() as workdir:
        workdir = Path(workdir)

        def run(name, options=(), engine=binary):
            proc, output = dock(engine, workdir, name, options)
            require(proc.returncode == 0, proc.stderr or proc.stdout)
            return proc.stdout, output.read_text()

        _, old = run("old", engine=baseline)
        _, default = run("default")
        _, native = run("native", ("--refinement_optimizer", "native"))
        require(old == default == native, "Default final refinement differs from the original optimizer")

        options = ("--refinement_optimizer", "armijo")
        _, armijo = run("armijo", options)
        require(len(parse(armijo)) == len(parse(old)), "Armijo changed the preselected refinement count")
        require(armijo != old, "The fixed fixture did not observe Armijo change the refined output")

        late = ("--pose_selection", "late")
        before_log, _ = run("native-late", late, engine=baseline)
        late_log, armijo_late = run("armijo-late", options + late)
        native_late_poses = parse(armijo_late)
        (before_counts,) = COUNTS.findall(before_log)
        (late_counts,) = COUNTS.findall(late_log)
        require(before_counts[:2] == late_counts[:2], "Armijo changed upstream late candidate counts")
        require(int(late_counts[3]) == len(native_late_poses), "Armijo late retained count differs from PDBQT")

        geometry = ("--refinement_optimizer", "geometry")
        geometry_log, geometry_late = run("geometry-late", geometry + late)
        (geometry_counts,) = COUNTS.findall(geometry_log)
        require(geometry_counts[:2] == late_counts[:2], "Geometry changed the upstream candidate counts")
        require(int(geometry_counts[3]) == len(parse(geometry_late)), "Geometry output count mismatch")

        pro_log, pro = run("armijo-late-pro", options + late + ("--rerank_gauss2_scale", "0.5"))
        pro_poses = sorted(parse(pro, scale=0.5), key=lambda pose: pose[2])
        require(COUNTS.findall(pro_log) == [late_counts], "Gauss2 changed Armijo refinement counts")
        native_records = [(p[0], p[2], p[3]) for p in native_late_poses]
        require([(p[0], p[2], p[3]) for p in pro_poses] == native_records,
                "Gauss2 changed Armijo-refined configurations or native scores")

        for name, flags, message in (
            ("invalid-mode", ("--refinement_optimizer", "other"), "must be native, armijo, geometry, or rigid_metric"),
            ("no-refine", options + ("--no_refine",), "Armijo refinement requires"),
            ("maps", options + ("--maps", "unused"), "Armijo refinement requires"),
            ("local-score-only", geometry + ("--local_only", "--score_only"), "Armijo refinement requires"),
            ("local-no-refine", options + ("--local_only", "--no_refine"), "Armijo refinement requires"),
            ("rigid-metric-global", ("--refinement_optimizer", "rigid_metric"), "rigid_metric requires local_only"),
            ("vinardo", options + ("--scoring", "vinardo"), "Armijo refinement requires"),
        ):
            proc, output = dock(binary, workdir, name, flags)
            require(proc.returncode != 0 and message in proc.stderr and not output.exists(),
                    f"Unexpected handling for {name}: {proc.stderr}")

        check_local_budget(binary, workdir)

    print("Native byte parity, Armijo/geometry, Late/Gauss2, local gradient caps and invalid modes passed.")


if __name__ == "__main__":
    main()
