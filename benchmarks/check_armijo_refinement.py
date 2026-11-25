#!/usr/bin/env python3
"""Check Armijo opt-in behavior, default parity and composition on a real Vina case."""

import argparse
from pathlib import Path
import tempfile

from check_late_pose_selection import COUNTS, dock, parse, require


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

        pro_log, pro = run("armijo-late-pro", options + late + ("--rerank_gauss2_scale", "0.5"))
        pro_poses = sorted(parse(pro, scale=0.5), key=lambda pose: pose[2])
        require(COUNTS.findall(pro_log) == [late_counts], "Gauss2 changed Armijo refinement counts")
        native_records = [(p[0], p[2], p[3]) for p in native_late_poses]
        require([(p[0], p[2], p[3]) for p in pro_poses] == native_records,
                "Gauss2 changed Armijo-refined configurations or native scores")

        for name, flags, message in (
            ("invalid-mode", ("--refinement_optimizer", "other"), "must be native or armijo"),
            ("no-refine", options + ("--no_refine",), "Armijo refinement requires regular"),
            ("maps", options + ("--maps", "unused"), "Armijo refinement requires regular"),
            ("local", options + ("--local_only",), "Armijo refinement requires regular"),
            ("vinardo", options + ("--scoring", "vinardo"), "Armijo refinement requires regular"),
        ):
            proc, output = dock(binary, workdir, name, flags)
            require(proc.returncode != 0 and message in proc.stderr and not output.exists(),
                    f"Unexpected handling for {name}: {proc.stderr}")

    print("Native output byte parity, Armijo refinement, Late/Gauss2 composition and invalid modes passed.")


if __name__ == "__main__":
    main()
