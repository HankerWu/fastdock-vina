#!/usr/bin/env python3
"""Check polar reranking through the normal Vina CLI and output window."""

import argparse
from pathlib import Path
import math
import subprocess
from tempfile import TemporaryDirectory

from check_late_pose_selection import COUNTS, ROOT, parse, require


LATE = ("--pose_selection", "late")
MODES = ("distance", "directional", "penalty")
HEADERS = ("REMARK VINA POLAR RESULT:", "REMARK VINA NATIVE:", "REMARK VINA POLAR SCORES:",
           "REMARK VINA POLAR TERMS:", "REMARK VINA POLAR MODE:")


def read_polar(text, mode, n_tors):
    poses = []
    for block in text.split("ENDMDL"):
        lines = block.strip().splitlines()
        if not lines:
            continue
        require(lines[0] == f"MODEL {len(poses) + 1}", "Output MODELs do not follow ranking")
        require(tuple(line.partition(":")[0] + ":" for line in lines[1:6]) == HEADERS,
                "Polar PDBQT headers differ from their contract")
        selected, distance, directional, radial, angular, native, rank, lower, upper = decode(lines, mode)
        g = 1 / (1 + n_tors * 0.05846)
        require(0 <= angular <= radial and math.isclose(distance, native - g * radial, rel_tol=0, abs_tol=1e-12)
                and math.isclose(directional, native - g * angular, rel_tol=0, abs_tol=1e-12),
                "Polar scores do not follow the original torsion transform")
        penalty = native + g * (radial - angular)
        require(math.isclose(penalty, native + directional - distance, rel_tol=0, abs_tol=1e-12),
                "Angular penalty is inconsistent with the two reported ranking scores")
        expected = {"distance": distance, "directional": directional, "penalty": penalty}
        require(math.isclose(selected, expected[mode], rel_tol=0, abs_tol=1e-12),
                "Selected rank score does not match the requested mode")
        diagnostic = (native, distance, directional, radial, angular)
        atoms = tuple(line for line in lines if line.startswith(("ATOM  ", "HETATM")))
        poses.append((selected, rank, atoms, diagnostic, lower, upper))
    require(poses and all(a[:2] <= b[:2] for a, b in zip(poses, poses[1:])), "Polar ranks are out of order")
    require(poses[0][-2:] == (0, 0), "RMSD must be relative to the first emitted pose")
    return poses


def decode(lines, mode):
    require(lines[5] == HEADERS[4] + " " + mode, "Wrong polar mode is written")
    require(sum(line.startswith((*HEADERS, "REMARK VINA RESULT:", "REMARK VINA PRO RESULT:"))
                for line in lines) == len(HEADERS), "Stale or misleading score header in a polar pose")
    selected, lower, upper = map(float, lines[1].split(":")[1].split())
    native, rank = lines[2].split(":")[1].split()
    distance, directional = map(float, lines[3].split(":")[1].split())
    radial, angular = map(float, lines[4].split(":")[1].split())
    require(all(map(math.isfinite, (selected, float(native), distance, directional, radial, angular))),
            "Nonfinite polar ranking diagnostic")
    return selected, distance, directional, radial, angular, float(native), int(rank), lower, upper


def match_native(guided, native):
    require(sorted(p[1] for p in guided) == list(range(1, len(native) + 1)), "Native pose ranks changed")
    for _, rank, atoms, diagnostics, *_ in guided:
        saved = native[rank - 1]
        require(atoms == saved[0] and abs(saved[1] - diagnostics[0]) <= 0.0005001,
                "Polar rank changed a native pose or score")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path, help="verified engine before introducing polar ranking")
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    baseline = args.baseline.resolve(strict=True)
    fixture = ROOT / "examples/basic"
    ligand = fixture / "1iep_ligand.pdbqt"
    with TemporaryDirectory(prefix="vina-polar-check-") as temp:
        directory = Path(temp)

        def dock(name, options=(), engine=binary, energy_range=100, input_ligand=ligand, error=None):
            output = directory / f"{name}.pdbqt"
            command = [str(engine), "--receptor", str(fixture / "1iep_receptor.pdbqt"),
                       "--ligand", str(input_ligand), "--config", str(fixture / "1iep_receptor.box.txt"),
                       "--out", str(output), "--cpu", "1", "--seed", "42", "--exhaustiveness", "2",
                       "--max_evals", "1500", "--num_modes", "5", "--min_rmsd", "1", "--energy_range", str(energy_range)]
            result = subprocess.run(command + list(options), capture_output=True, text=True, timeout=90)
            if error:
                require(result.returncode != 0 and error in result.stderr and not output.exists(),
                        f"{name}: unexpected result {result.stderr}")
                return None
            require(result.returncode == 0, result.stderr or result.stdout)
            return result.stdout, output.read_text()

        for flags, name in (((), "native"), (LATE, "late")):
            old, previous = dock(name + "-old", flags, engine=baseline)
            log, current = dock(name, flags)
            require(current == previous and COUNTS.findall(log) == COUNTS.findall(old),
                    f"Default {name} changed with polar ranking disabled")
        native = parse(current)
        native_counts = COUNTS.findall(log)
        full = {}
        for mode in MODES:
            options = LATE + ("--polar_geometry", mode)
            log, output = dock(mode, options)
            require(COUNTS.findall(log) == native_counts, "Polar rank changed preceding Late stages")
            guide = read_polar(output, mode, n_tors=7)
            match_native(guide, native)
            require(any(p[3][-1] > 0 for p in guide), "Fixture did not exercise angular scoring")
            full[mode] = {p[1]: p[2:4] for p in guide}
            require(f"Polar geometry scoring: mode={mode} poses=" in log, "Missing polar work count")
            _, tight_text = dock(mode + "-tight", options, energy_range=0)
            tight = read_polar(tight_text, mode, n_tors=7)
            best = min(p[3][0] for p in guide)
            require([p[:4] for p in tight] == [p[:4] for p in guide if p[3][0] == best],
                    "Native energy window changed after polar reranking")
        require(full["distance"] == full["directional"] == full["penalty"],
                "Polar modes disagree on common candidates or raw diagnostics")

        annotated = directory / "annotated.pdbqt"
        old_headers = "\n".join(output.splitlines()[1:6])
        annotated.write_text(old_headers + "\nREMARK TEST NON-SCORING NOTE\n" + ligand.read_text())
        _, cleaned = dock("cleanup", LATE + ("--polar_geometry", "directional"), input_ligand=annotated)
        require("REMARK TEST NON-SCORING NOTE" in cleaned, "A non-score input remark was discarded")
        read_polar(cleaned, "directional", n_tors=7)

        carbon = next(line for line in ligand.read_text().splitlines() if line.startswith("ATOM  ")
                      and line[77:].strip() == "C")
        one_atom = directory / "unpolar.pdbqt"
        one_atom.write_text(f"ROOT\n{carbon}\nENDROOT\nTORSDOF 0\n")
        _, plain = dock("zero-native", LATE, input_ligand=one_atom)
        _, inactive = dock("zero-polar", LATE + ("--polar_geometry", "directional"), input_ligand=one_atom)
        guided = read_polar(inactive, "directional", n_tors=0)
        match_native(guided, parse(plain))
        require(all(p[1] == i and p[3][-2:] == (0, 0) and p[0] == p[3][0]
                    for i, p in enumerate(guided, 1)), "No-pair pool changed score or order")

        for name, options, message in (
            ("not-late", ("--polar_geometry", "distance"), "polar_geometry requires"),
            ("bad-mode", LATE + ("--polar_geometry", "invalid"), "must be distance, directional or penalty"),
            ("gauss2", LATE + ("--polar_geometry", "directional", "--rerank_gauss2_scale", "0.5"), "polar_geometry requires"),
            ("template", LATE + ("--polar_geometry", "penalty", "--pose_template", "unused", "--template_stage", "rank"),
             "polar_geometry requires"),
        ):
            dock(name, options, error=message)
    print("Polar CLI: native parity, same pools, three scores, zero-pair fallback and native energy window passed.")


if __name__ == "__main__":
    main()
