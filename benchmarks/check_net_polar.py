#!/usr/bin/env python3
"""Validate the net polar scorer's analytic terms and native-pool parity."""

import argparse
from math import pi, sqrt
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from check_late_pose_selection import COUNTS, ROOT, dock, parse, require
from check_polar_burial import atom, same, values


def fixed(binary, receptor, ligand, options=(), mode="combined", check=True):
    command = [str(binary), "--receptor", str(receptor), "--ligand", str(ligand),
               "--score_only", "--net_polar", mode, "--cpu", "1", "--verbosity", "0"]
    return subprocess.run(command + list(options), capture_output=True, text=True, timeout=30, check=check)


def diagnostics(text):
    terms, deltas, native, (normalization,), (water_points,) = (
        values(text, f"Net polar {label}:")
        for label in ("terms", "deltas", "native", "normalization", "water points")
    )
    _, _, direction, buried_area, water_area = terms
    raw_burial = 0.02 * buried_area
    expected = (direction, raw_burial, direction + raw_burial,
                direction + raw_burial - 0.02 * water_area)
    require(same(deltas, tuple(value * normalization for value in expected)), "Wrong reported score deltas")
    require(same((native[2],), ((native[0] + native[1]) * normalization,)),
            "Wrong common-reference normalization")
    require(0 <= water_area <= buried_area and (water_points > 0 or water_area == 0),
            "Water credit outside the permitted burial area")
    return terms, deltas, native, normalization


def check_analytic(binary, directory):
    receptor, ligand, waters = [directory / name for name in ("receptor.pdbqt", "ligand.pdbqt", "waters.txt")]

    def write_ligand(atoms):
        ligand.write_text("ROOT\n" + atoms + "ENDROOT\nTORSDOF 0\n")

    # Equal expanded oxygen radii at one radius distance occlude one quarter
    # of the ligand sphere. Fibonacci z coordinates hit this cap exactly.
    write_ligand(atom(1, "OA", 0))
    receptor.write_text(atom(1, "OA", 0, z=2.92))
    terms, _, _, normalization = diagnostics(fixed(binary, receptor, ligand).stdout)
    require(same(terms, (0, 0, 0, pi * 2.92**2, 0)) and normalization == 1, "Wrong analytic oxygen cap")
    # This covalent carbon already occludes the entire receptor cap.
    write_ligand(atom(1, "OA", 0) + atom(2, "C", 0, z=1.3))
    require(same(diagnostics(fixed(binary, receptor, ligand).stdout)[0], (0, 0, 0, 0, 0)),
            "Area already occluded in the isolated ligand was counted")
    # A distant receptor contributes neither direction nor area.
    receptor.write_text(atom(1, "OA", 0, z=9))
    require(same(diagnostics(fixed(binary, receptor, ligand).stdout)[0], (0, 0, 0, 0, 0)),
            "Noninteracting receptor received a polar term")

    # Native term uses the exact midpoint endpoints in r², with a custom HB
    # weight to distinguish reuse of the native term from an approximate one.
    write_ligand(atom(1, "OA", 0) + atom(2, "C", 1.3))
    distance, weight = 3.15, -1.13
    bin_index = int(32 * distance**2)
    native_hb = weight * (3.5 - (sqrt(bin_index / 32) + sqrt((bin_index + 1) / 32)) / 2) / 0.7
    for sign in (-1, +1):
        receptor.write_text(atom(1, "N", sign * distance) + atom(2, "HD", sign * distance, z=1.02))
        terms = diagnostics(fixed(binary, receptor, ligand, ("--weight_hydrogen", str(weight))).stdout)[0]
        directed = native_hb if sign < 0 else 0
        require(same(terms[:2], (native_hb, directed)), "Wrong native or replaced hydrogen bond")
        require(terms[2] == 0 if sign < 0 else 0 < terms[2] <= -native_hb,
                "Wrong replacement after original non-linear scoring")

    # Water is aligned at (2.9, 0, 0), giving two ideal distances and a 90°
    # bridge while the receptor still occludes some ligand oxygen surface.
    write_ligand(atom(1, "OA", 0) + atom(2, "C", -1.3))
    receptor.write_text(atom(1, "OA", 2.9, y=2.9) + atom(2, "C", 2.9, y=4.2))
    no_water = fixed(binary, receptor, ligand, mode="water", check=False)
    require(no_water.returncode != 0 and "independent aligned structural evidence" in no_water.stderr,
            "Water arm accepted missing evidence")
    waters.write_text("# independent aligned water O\n2.9 0 0 # one oxygen\n")
    flags = ("--net_polar_waters", str(waters))
    terms, deltas, _, _ = diagnostics(fixed(binary, receptor, ligand, flags, mode="water").stdout)
    require(terms[3] > 0 and same((terms[4], deltas[3]), (terms[3], 0)), "Ideal bridge did not offset the area penalty")
    waters.write_text("0 0 0\n")
    terms, deltas, _, _ = diagnostics(fixed(binary, receptor, ligand, flags, mode="water").stdout)
    require(terms[3] > 0 and terms[4] == 0 and deltas[2] == deltas[3], "Clashing water gave compensation")
    waters.write_text("2.9 0 0 extra\n")
    malformed = fixed(binary, receptor, ligand, flags, mode="water", check=False)
    require(malformed.returncode != 0 and "three finite coordinates" in malformed.stderr,
            "Malformed water evidence was accepted")


def check_real(binary, baseline, directory):
    late = ("--pose_selection", "late", "--refine_candidates", "6")
    old, previous = dock(baseline, directory, "old", late)
    new, current = dock(binary, directory, "new", late)
    old.check_returncode()
    new.check_returncode()
    require(previous.read_bytes() == current.read_bytes() and
            COUNTS.findall(old.stdout) == COUNTS.findall(new.stdout), "Default Late output changed")
    native = parse(current.read_text())
    fixture = ROOT / "examples/basic"
    receptor = fixture / "1iep_receptor.pdbqt"
    ligand = fixture / "1iep_ligand.pdbqt"
    report = fixed(binary, receptor, ligand).stdout
    terms, _, native_fixed, scale = diagnostics(report)
    require(terms[3] > 0 and abs(scale - 1 / (1 + 7 * 0.05846)) < 1e-12,
            "Real input did not use Vina's actual normalization")
    previous_report = subprocess.run(
        [str(baseline), "--receptor", str(receptor), "--ligand", str(ligand), "--score_only",
         "--polar_burial", "compensated", "--cpu", "1", "--verbosity", "0"],
        check=True, capture_output=True, text=True, timeout=30,
    ).stdout
    require(native_fixed == values(previous_report, "Polar burial native:"), "Native fixed score changed")

    common = None
    for mode in ("native", "direction", "burial", "combined"):
        run, output = dock(binary, directory, mode, late + ("--net_polar", mode))
        run.check_returncode()
        require(COUNTS.findall(run.stdout) == COUNTS.findall(new.stdout), "Search or refinement budget changed")
        order, raw_by_rank = [], {}
        for block in filter(str.strip, output.read_text().split("ENDMDL")):
            energy, rank_float = values(block, "REMARK VINA NATIVE:")
            rank = int(rank_float)
            selected, lower, upper = values(block, "REMARK VINA NET POLAR RESULT:")
            terms = values(block, "REMARK VINA NET POLAR TERMS:")
            deltas = values(block, "REMARK VINA NET POLAR DELTAS:")
            atoms = tuple(line for line in block.splitlines() if line.startswith(("ATOM  ", "HETATM")))
            require(atoms == native[rank - 1][0], f"Native pose {rank} changed coordinates")
            require(abs(energy - native[rank - 1][1]) < 0.000501, "Native affinity changed")
            expected = energy if mode == "native" else energy + deltas[("direction", "burial", "combined").index(mode)]
            require(same((selected,), (expected,)), "Selected full-precision score does not match the requested arm")
            if not order:
                require((lower, upper) == (0, 0), "RMSD reference does not match the first emitted pose")
            order.append((selected, rank))
            raw_by_rank[rank] = (terms, deltas)
        require(order == sorted(order) and sorted(raw_by_rank) == list(range(1, len(native) + 1)),
                "Ranking is not sorted or retained a different pool")
        if common is None:
            common = raw_by_rank
        else:
            require(raw_by_rank == common, "Raw score differs across arms on the same pose")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    baseline = args.baseline.resolve(strict=True)
    with TemporaryDirectory(prefix="vina-net-polar-check-") as temp:
        directory = Path(temp)
        check_analytic(binary, directory)
        check_real(binary, baseline, directory)
    print("Net polar analytic SASA, native H-bond replacement, evidence water, and 1IEP native-pool parity passed.")


if __name__ == "__main__":
    main()
