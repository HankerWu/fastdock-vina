#!/usr/bin/env python3
"""Check the fixed-coordinate scorer and native-pool preservation."""

import argparse
from math import isclose
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from check_late_pose_selection import COUNTS, ROOT, dock, parse, require


MODES = ("burial", "geometry", "compensated")
LATE = ("--pose_selection", "late", "--refine_candidates", "6")


def values(text, prefix):
    (line,) = [line for line in text.splitlines() if line.startswith(prefix)]
    return tuple(map(float, line[len(prefix):].split()))


def same(first, second, tolerance=1e-11):
    return len(first) == len(second) and all(
        isclose(x, y, abs_tol=tolerance, rel_tol=0) for x, y in zip(first, second)
    )


def fixed(binary, rec, lig):
    process = subprocess.run(
        [str(binary), "--receptor", str(rec), "--ligand", str(lig), "--score_only",
         "--polar_burial", "compensated", "--cpu", "1", "--verbosity", "0"],
        capture_output=True, text=True, timeout=30, check=True,
    )
    return tuple(values(process.stdout, "Polar burial " + key + ":")
                 for key in ("terms", "deltas", "native"))


def atom(serial, kind, x, y=0, z=0):
    return (f"ATOM  {serial:5d}  X   UNL     1    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00     0.000 {kind}\n")


def check_analytic(binary, directory):
    receptor, ligand = directory / "rigid.pdbqt", directory / "ligand.pdbqt"
    ligand.write_text("ROOT\n" + atom(1, "OA", 0) + atom(2, "C", 1.30) + "ENDROOT\nTORSDOF 0\n")
    cases = (
        (atom(1, "C", 5.5), (0, 0, 0)),
        (atom(1, "C", -4.5), (1 / 16, 0, 1 / 16)),
        (atom(1, "C", -2.9), (1 / 8, 0, 1 / 8)),
        (atom(1, "Zn", -2.9), (1 / 8, 0, 1 / 8)),
        (atom(1, "N", -2.9) + atom(2, "HD", -2.9, z=1.02), (1 / 8, 0, 0)),
        (atom(1, "N", 2.9) + atom(2, "HD", 2.9, z=1.02), (1 / 8, 1, 1 / 8)),
        ("".join(atom(i + 1, "C", -2, y=i / 4) for i in range(9)), (1, 0, 1)),
    )
    for atoms, expected in cases:
        receptor.write_text(atoms)
        terms, deltas, native = fixed(binary, receptor, ligand)
        require(same(terms, expected), f"Expected raw terms {expected}, got {terms}")
        require(same(deltas, tuple(0.6 * value for value in terms)), "Incorrect zero-torsion penalty")
        require(same((native[0] + native[1],), (native[2],)), "Incorrect diagnostic common reference")


def check_pool(binary, baseline, directory):
    old, previous = dock(baseline, directory, "old", LATE)
    new, current = dock(binary, directory, "new", LATE)
    old.check_returncode()
    new.check_returncode()
    require(previous.read_bytes() == current.read_bytes() and
            COUNTS.findall(old.stdout) == COUNTS.findall(new.stdout), "Default Late changed")
    native = parse(current.read_text())
    receptor = ROOT / "examples/basic/1iep_receptor.pdbqt"
    pose = directory / "fixed.pdbqt"
    all_terms = {}
    for arm_index, mode in enumerate(MODES):
        run, output = dock(binary, directory, mode, LATE + ("--polar_burial", mode))
        run.check_returncode()
        require(COUNTS.findall(run.stdout) == COUNTS.findall(new.stdout), "Search or refinement budget changed")
        order, treatment = [], {}
        for model in filter(str.strip, output.read_text().split("ENDMDL")):
            native_score, raw_rank = values(model, "REMARK VINA NATIVE:")
            rank = int(raw_rank)
            rank_score, lower, upper = values(model, "REMARK VINA POLAR BURIAL RESULT:")
            terms = values(model, "REMARK VINA POLAR BURIAL TERMS:")
            deltas = values(model, "REMARK VINA POLAR BURIAL DELTAS:")
            atoms = tuple(line for line in model.splitlines() if line.startswith(("ATOM  ", "HETATM")))
            require(atoms == native[rank - 1][0], f"Pose {rank} changed coordinates")
            require(abs(native_score - native[rank - 1][1]) < 0.000501, "Native affinity changed")
            require(same(deltas, tuple(0.6 * value / (1 + 7 * 0.05846) for value in terms)),
                    "Incorrect Vina torsion normalization")
            require(same((rank_score,), (native_score + deltas[arm_index],)), "Incorrect rank score")
            require(0 <= terms[2] <= terms[0] and terms[1] >= 0, "Invalid raw terms")
            if not order:
                require((lower, upper) == (0, 0), "Reference RMSD does not match first emitted pose")
            order.append((rank_score, rank))
            treatment[rank] = (terms, deltas)
            if arm_index == 0:
                pose.write_text("\n".join(line for line in model.splitlines()
                                           if line and not line.startswith("MODEL")) + "\n")
                fixed_terms, fixed_deltas, _ = fixed(binary, receptor, pose)
                require(same(fixed_terms, terms, 0.025) and same(fixed_deltas, deltas, 0.011),
                        "Fixed coordinates disagree beyond expected PDBQT rounding")
        require(order == sorted(order) and sorted(treatment) == list(range(1, len(native) + 1)),
                "Ranking is not sorted or retained a different pool")
        if arm_index == 0:
            all_terms = treatment
        else:
            require(all_terms == treatment, "Arms disagree on shared raw terms")
    rejected, output = dock(binary, directory, "rejected", ("--polar_burial", "burial"))
    require(rejected.returncode != 0 and not output.exists(), "Accepted non-Late search")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    baseline = args.baseline.resolve(strict=True)
    with TemporaryDirectory(prefix="vina-polar-burial-") as temp:
        directory = Path(temp)
        check_analytic(binary, directory)
        check_pool(binary, baseline, directory)
    print("Polar burial numerical cases, stable native pool, and fixed-coordinate parity passed.")


if __name__ == "__main__":
    main()
