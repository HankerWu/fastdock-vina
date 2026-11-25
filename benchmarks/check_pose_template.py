#!/usr/bin/env python3
"""Check native-mode parity, template ranking and the three docking entrypoints."""

import argparse
import math
from pathlib import Path
import re
import subprocess
import tempfile
import time

from check_late_pose_selection import COUNTS, ROOT, parse, require


FIXTURE = ROOT / "examples/basic"
TEMPLATE_COUNTS = re.compile(
    r"Template pose scoring: stage=(\w+) selection_evals=(\d+) refinement_evals=(\d+) rank_evals=(\d+)"
)
RESULT = "REMARK VINA TEMPLATE RESULT:"
NATIVE = "REMARK VINA NATIVE:"
OVERLAP = "REMARK VINA TEMPLATE OVERLAP:"
STAGE = "REMARK VINA TEMPLATE STAGE:"
LATE = ("--pose_selection", "late")
ELEMENTS = {"A": "C", "NA": "N", "OA": "O", "SA": "S", **{f"CG{i}": "C" for i in range(4)}}
EXCLUDED = {"H", "HD", "W", *(f"G{i}" for i in range(4))}


def real_atoms(lines):
    atoms = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            kind = line[77:].strip()
            if kind not in EXCLUDED:
                atoms.append((ELEMENTS.get(kind, kind),
                              tuple(float(line[i:i + 8]) for i in (30, 38, 46))))
    return atoms


def gaussian_overlap(ligand, template):
    def directional(first, second):
        return sum(max((math.exp(-math.dist(a, b) ** 2 / 4.5)
                        for matched, b in second if element == matched), default=0)
                   for element, a in first) / len(first)

    return (directional(ligand, template) + directional(template, ligand)) / 2


def read_guided(text, stage, template):
    poses = []
    for block in text.split("ENDMDL"):
        if not block.strip():
            continue
        lines = block.strip().splitlines()
        require(lines[0] == f"MODEL {len(poses) + 1}", "Model numbers must follow template score")
        require(tuple(line.partition(":")[0] + ":" for line in lines[1:5]) == (RESULT, NATIVE, OVERLAP, STAGE),
                "Unexpected template score header")
        require(lines[4] == f"{STAGE} {stage}", "The template stage was not retained exactly")
        require(sum(line.startswith((RESULT, NATIVE, OVERLAP, STAGE, "REMARK VINA RESULT:",
                                     "REMARK VINA PRO RESULT:", "REMARK VINA PRO GAUSS2 SCALE:"))
                    for line in lines) == 4, "Output retains a stale score header")
        combined, low, high = map(float, lines[1].removeprefix(RESULT).split())
        native, native_rank = lines[2].removeprefix(NATIVE).split()
        native = float(native)
        overlap = float(lines[3].removeprefix(OVERLAP))
        atoms = tuple(line for line in lines if line.startswith(("ATOM  ", "HETATM")))
        require(0 <= overlap <= 1 and abs(overlap - gaussian_overlap(real_atoms(atoms), template)) < 0.0004,
                "PDBQT overlap cannot be independently reproduced from rounded output coordinates")
        # The example has seven qualifying Vina rotors at the standard torsion weight.
        require(abs(combined - (native - 6 * overlap / (1 + 7 * 0.05846))) < 1e-12,
                "Wrong Vina torsion transformation for the reported ranking score")
        poses.append((atoms, combined, int(native_rank), native, low, high))
    require(poses and all(a[1:3] <= b[1:3] for a, b in zip(poses, poses[1:])), "Template poses are not sorted")
    require(poses[0][-2:] == (0, 0), "Distances must be relative to the first emitted pose")
    return poses


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path, help="verified build before template support")
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    baseline = args.baseline.resolve(strict=True)
    ligand = FIXTURE / "1iep_ligand.pdbqt"
    with tempfile.TemporaryDirectory(prefix="vina-template-check-") as directory:
        directory = Path(directory)

        def dock(name, options=(), engine=binary, energy_range=100, input_ligand=ligand, error=None):
            output = directory / f"{name}.pdbqt"
            command = [str(engine), "--receptor", str(FIXTURE / "1iep_receptor.pdbqt"),
                       "--ligand", str(input_ligand), "--config", str(FIXTURE / "1iep_receptor.box.txt"),
                       "--out", str(output), "--cpu", "1", "--exhaustiveness", "3", "--max_evals", "1800",
                       "--num_modes", "3", "--min_rmsd", "1", "--energy_range", str(energy_range),
                       "--seed", "42", *options]
            start = time.perf_counter()
            completed = subprocess.run(command, capture_output=True, text=True, timeout=90)
            elapsed = time.perf_counter() - start
            if error:
                require(completed.returncode != 0 and error in completed.stderr and not output.exists(),
                        f"{name}: unexpected result: {completed.stderr}")
                return
            require(completed.returncode == 0, completed.stderr or completed.stdout)
            return completed.stdout, output.read_text(), elapsed

        for name, flags in (("default", ()), ("late", LATE),
                            ("armijo", LATE + ("--refinement_optimizer", "armijo")),
                            ("gauss2", LATE + ("--rerank_gauss2_scale", "0.5"))):
            old_log, old_result, _ = dock(name + "-old", flags, engine=baseline)
            log, result, elapsed = dock(name, flags)
            require(result == old_result and COUNTS.findall(log) == COUNTS.findall(old_log),
                    f"Existing {name} mode changed from the verified build")
            if name == "late":
                native_text = result
                native_counts = tuple(map(int, COUNTS.findall(log)[0]))
                native_time = elapsed

        native_poses = parse(native_text)
        # A known pool member is used only to make this engineering check sensitive to reranking.
        template_atoms = real_atoms(native_poses[-1][0])
        guide = directory / "template.txt"
        guide.write_text("".join(f"{e} {x} {y} {z}\n" for e, (x, y, z) in template_atoms))
        template_flags = ("--pose_template", str(guide))
        for stage in ("rank", "select", "refine"):
            flags = LATE + template_flags + ("--template_stage", stage)
            log, output, elapsed = dock(stage, flags)
            poses = read_guided(output, stage, template_atoms)
            ((saved, inputs, calls, retained),) = [tuple(map(int, m)) for m in COUNTS.findall(log)]
            ((mode, selected_evals, refined_evals, ranked_evals),) = TEMPLATE_COUNTS.findall(log)
            require("template score" in log and mode == stage and len(poses) == retained,
                    "Missing diagnostics about the scored candidate pool")
            require(saved == native_counts[0] and inputs == native_counts[1] and inputs <= calls <= inputs * 5,
                    "Changed trajectory budget or refinement cap")
            require(int(selected_evals) == (0 if stage == "rank" else saved), "Incorrect selection work count")
            require((int(refined_evals) > 0) == (stage == "refine"), "Incorrect BFGS template work count")
            require(int(ranked_evals) == (inputs if stage == "refine" else retained), "Incorrect final overlap work count")
            if stage == "rank":
                require(tuple(row[0] for row in sorted(poses, key=lambda row: row[2]))
                        == tuple(p[0] for p in native_poses), "Rank mode changed candidate coordinates or native order")
                require(all(abs(row[3] - native_poses[row[2] - 1][1]) <= 0.0005 for row in poses),
                        "Rank mode changed native scores")
                require([row[2] for row in poses] == [1, 3, 2], "The fixture did not exercise template-driven reorder")
                window = (native_poses[1][1] + native_poses[2][1]) / 2 - native_poses[0][1]
                _, tight_text, _ = dock("rank-tight", flags, energy_range=window)
                tight = read_guided(tight_text, stage, template_atoms)
                _, native_tight, _ = dock("late-tight", LATE, energy_range=window)
                require(len(tight) == 2 and {row[0] for row in tight} == {row[0] for row in parse(native_tight)},
                        "Rank mode changed native energy-window membership")
            print(f"{stage}: {elapsed:.3f} s, Late {saved}/{inputs}/{calls}/{retained}, "
                  f"template evaluations {selected_evals}/{refined_evals}/{ranked_evals}")
        print(f"Native Late fixture: {native_time:.3f} s")

        stale = directory / "stale.pdbqt"
        stale.write_text("\n".join((f"{RESULT} 0 0 0", f"{NATIVE} 0 1", f"{OVERLAP} 0",
                                    f"{STAGE} select", "REMARK VINA PRO GAUSS2 SCALE: 0.5",
                                    "REMARK TEST NON-SCORING NOTE")) + "\n" + ligand.read_text())
        _, clean, _ = dock("rank-stale", LATE + template_flags + ("--template_stage", "rank"), input_ligand=stale)
        read_guided(clean, "rank", template_atoms)
        require("REMARK TEST NON-SCORING NOTE" in clean, "A non-score remark was lost")

        invalid = directory / "invalid.txt"
        invalid.write_text("C 0 1 broken\n")
        empty = directory / "empty.txt"
        empty.write_text("H 0 0 0\nG0 1 2 3\n")
        for name, flags, error in (
            ("unpaired-template", LATE + template_flags, "Specify pose_template and template_stage together"),
            ("unpaired-stage", LATE + ("--template_stage", "rank"), "Specify pose_template and template_stage together"),
            ("bad-stage", LATE + template_flags + ("--template_stage", "invalid"), "must be rank, select or refine"),
            ("without-late", template_flags + ("--template_stage", "rank"), "pose_template requires Late Vina"),
            ("armijo-template", LATE + template_flags + ("--template_stage", "rank", "--refinement_optimizer", "armijo"),
             "pose_template requires Late Vina"),
            ("gauss-template", LATE + template_flags + ("--template_stage", "rank", "--rerank_gauss2_scale", "0.5"),
             "pose_template requires Late Vina"),
            ("bad-coordinate", LATE + ("--pose_template", str(invalid), "--template_stage", "refine"),
             "Expected element and three finite coordinates"),
            ("no-heavy-atoms", LATE + ("--pose_template", str(empty), "--template_stage", "rank"),
             "must each contain real heavy atoms"),
        ):
            dock(name, flags, error=error)


if __name__ == "__main__":
    main()
