#!/usr/bin/env python3
"""Verify local explicit refinement and retention of grid modes."""

import argparse
from pathlib import Path
import re
import subprocess
import tempfile

from check_late_pose_selection import ROOT, require


AFFINITY = re.compile(r"Estimated Free Energy of Binding\s+:\s+([-\d.]+)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    baseline = args.baseline.resolve(strict=True)
    fixture = ROOT / "examples/basic"

    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)

        def call(engine, sf, ligand, options):
            command = [str(engine), "--scoring", sf, "--ligand", str(ligand), "--seed", "42", "--cpu", "1"]
            result = subprocess.run(command + options, text=True, capture_output=True, timeout=90)
            require(result.returncode == 0, result.stderr or result.stdout)
            (score,) = AFFINITY.findall(result.stdout)
            return float(score)

        def optimize(engine, sf, ligand, label, options):
            output = directory / f"{label}.pdbqt"
            score = call(engine, sf, ligand, options + ["--local_only", "--out", str(output)])
            atoms = tuple(line for line in output.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM")))
            require(atoms, "No atoms in the optimized ligand")
            return score, atoms, output

        receptor = ["--receptor", str(fixture / "1iep_receptor.pdbqt"),
                    "--config", str(fixture / "1iep_receptor.box.txt")]
        ligand = fixture / "1iep_ligand.pdbqt"
        for sf in ("vina", "vinardo"):
            old = optimize(baseline, sf, ligand, f"{sf}-old", receptor)
            old_grid = optimize(baseline, sf, ligand, f"{sf}-old-grid", receptor + ["--no_refine"])
            fixed = optimize(binary, sf, ligand, f"{sf}-explicit", receptor)
            maps = directory / f"{sf}-maps"
            fixed_grid = optimize(binary, sf, ligand, f"{sf}-grid",
                                  receptor + ["--no_refine", "--force_even_voxels", "--write_maps", str(maps)])
            require(old[1] == old_grid[1] and old[0] != old_grid[0],
                    "Baseline should reproduce grid optimization with explicit scoring")
            require(fixed[1] != old[1] and fixed[0] < old[0],
                    "Explicit receptor optimization did not change the 1IEP result")
            require(fixed_grid[:2] == old_grid[:2], "Direct receptor grid behavior changed")
            rescored = call(binary, sf, fixed[2], receptor + ["--score_only"])
            # PDBQT rounds coordinates to 0.001 Å before the second process reads them.
            require(abs(rescored - fixed[0]) <= 0.02, "Explicit optimization and rescoring disagree")
            for engine, label in ((baseline, "old"), (binary, "fixed")):
                result = optimize(engine, sf, ligand, f"{sf}-map-only-{label}", ["--maps", str(maps)])
                if label == "old":
                    original_map_result = result[:2]
                else:
                    require(result[:2] == original_map_result, "Loading maps changed local refinement")
            print(f"{sf} explicit affinity {old[0]:.3f} -> {fixed[0]:.3f}; grid {fixed_grid[0]:.3f}.")

        ligand = directory / "one-carbon.pdbqt"
        ligand.write_text("ROOT\nATOM      1  X   UNL     1       1.000   0.000   0.000  1.00  0.00"
                          "     0.000 C\nENDROOT\nTORSDOF 0\n")
        maps = directory / "ad4-maps"
        header = "GRID_PARAMETER_FILE fixture.gpf\nGRID_DATA_FILE fixture.fld\nMACROMOLECULE fixture.pdbqt\n"
        header += "SPACING 1.0\nNELEMENTS 4 4 4\nCENTER 0.0 0.0 0.0\n"
        for suffix in ("C", "e", "d"):
            values = [f"{x * x if suffix == 'C' else 0}\n" for _ in range(25) for x in range(-2, 3)]
            maps.with_suffix(f".{suffix}.map").write_text(header + "".join(values))
        args_ad4 = ["--maps", str(maps)]
        old_ad4 = optimize(baseline, "ad4", ligand, "ad4-old", args_ad4)
        fixed_ad4 = optimize(binary, "ad4", ligand, "ad4-fixed", args_ad4)
        require(old_ad4[:2] == fixed_ad4[:2], "AD4 map optimization or scoring changed")
        ligand.write_text(ligand.read_text().replace("   1.000", " 100.000"))
        for op in ("--score_only", "--local_only"):
            result = subprocess.run([binary, "--scoring", "ad4", "--ligand", ligand, *args_ad4,
                                     op, "--out", directory / "outside.pdbqt"], capture_output=True, text=True)
            require(result.returncode != 0 and "outside the grid box" in result.stderr + result.stdout,
                    "AD4 fixed-pose and local entry points must check the loaded AD4 lattice")

    print("Explicit and grid optimization, exact map-only parity and AD4 parity passed.")


if __name__ == "__main__":
    main()
