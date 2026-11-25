#!/usr/bin/env python3
"""Check the new net-polar arms through the Vina and fixed-pool CLIs."""

import argparse
import json
from math import cos, exp, radians
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from check_late_pose_selection import COUNTS, ROOT, dock, require
from check_net_polar import diagnostics, fixed
from check_polar_burial import atom, same, values


SATISFACTION_REMARK = "REMARK VINA NET POLAR SATISFACTION:"


def check_report(text, prefix="Net polar"):
    terms = values(text, f"{prefix} terms:")
    deltas = values(text, f"{prefix} deltas:")
    direct, water, direct_cost, water_cost, dry_delta, wet_delta = values(text, f"{prefix} satisfaction:")
    area, legacy_water = terms[3:]
    require(0 <= direct <= area and 0 <= water <= legacy_water and direct + water <= area + 1e-11,
            "A reported credit exceeds its eligible area")
    require(deltas[0] <= wet_delta + 1e-11 and wet_delta <= dry_delta + 1e-11,
            "Water created a bonus or increased the residual penalty")
    require(same((dry_delta, wet_delta),
                 (deltas[2] - direct_cost, deltas[2] - direct_cost - water_cost)),
            "Score differs from the reported base and credits")
    return terms, (direct, water), (dry_delta, wet_delta)


def check_analytic(binary, directory):
    receptor, ligand, waters = [directory / name for name in ("receptor.pdbqt", "ligand.pdbqt", "waters.xyz")]
    flags = ("--net_polar_waters", str(waters))
    missing = fixed(binary, ROOT / "examples/basic/1iep_receptor.pdbqt",
                    ROOT / "examples/basic/1iep_ligand.pdbqt", mode="satisfied_water", check=False)
    require(missing.returncode != 0 and "independent aligned structural evidence" in missing.stderr,
            "Water arm accepted missing independent coordinates")

    def write_ligand(atoms):
        ligand.write_text("ROOT\n" + atoms + "ENDROOT\nTORSDOF 0\n")

    def donor(serial, x, y=0, z=0):
        return atom(serial, "N", x, y, z) + atom(serial + 1, "HD", x, y, z + 1.02)

    write_ligand(atom(1, "OA", 0) + atom(2, "C", -1.3))
    good_radial = exp(-0.5 * ((3.15 - 2.9) / 0.45) ** 2)
    for protein, quality in (
        (donor(1, 2.9), 1),
        (donor(1, 3.15), good_radial),
        (donor(1, -0.75, 2.80), 0),
        (donor(1, 4.001), 0),
        (atom(1, "OA", 2.9), 0),
    ):
        receptor.write_text(protein)
        terms, credits, deltas = check_report(fixed(binary, receptor, ligand, mode="satisfied").stdout)
        require(terms[3] > 0 and same(credits, (quality * terms[3], 0)) and deltas[0] == deltas[1],
                f"Wrong direct quality {quality} or dry fallback for {protein.strip()}")
    write_ligand(atom(1, "OA", 0))
    receptor.write_text(donor(1, 2.9))
    terms, credits, _ = check_report(fixed(binary, receptor, ligand, mode="satisfied").stdout)
    require(terms[3] > 0 and credits == (0, 0), "Missing acceptor geometry received credit")

    # Both candidates can pay the same acceptor's desolvation cost. At 90°
    # the direct leg is partial and the two optimal water legs beat it.
    write_ligand(atom(1, "OA", 0) + atom(2, "C", -1.3))
    waters.write_text("2.9 0 0\n2.9 0 0\n")
    water_partner = atom(3, "OA", 2.9, y=2.9) + atom(4, "C", 2.9, y=4.2)
    receptor.write_text(donor(1, 0, -2.9) + water_partner)
    report = fixed(binary, receptor, ligand, flags, mode="satisfied_water").stdout
    terms, credits, _ = check_report(report)
    q = cos(radians(80)) / (cos(radians(80)) - cos(radians(110)))
    require(terms[3] > 0 and same((terms[4], *credits), (terms[3], q * terms[3], (1 - q) * terms[3])),
            "Water and direct qualities did not compete at the same acceptor")
    # Change the direction so the direct leg is optimal. The same water can
    # still bridge, but now it contributes zero additional credit.
    write_ligand(atom(1, "OA", 0) + atom(2, "C", -0.919, 0.919))
    terms, credits, _ = check_report(fixed(binary, receptor, ligand, flags, mode="satisfied_water").stdout)
    require(terms[3] > 0 and same((terms[4], *credits), (terms[3], terms[3], 0)),
            "An already satisfied site received a second credit")
    waters.write_text("0 0 0\n")
    terms, credits, deltas = check_report(fixed(binary, receptor, ligand, flags, mode="satisfied_water").stdout)
    require(terms[4] == credits[1] == 0 and deltas[0] == deltas[1],
            "A clashing water changed the dry score")

    # Water pays for an acceptor whose direct credit is zero, even when a
    # second acceptor elsewhere is fully satisfied by direct protein contact.
    write_ligand(atom(1, "OA", 0) + atom(2, "C", -1.3)
                 + atom(3, "OA", 20) + atom(4, "C", 18.7))
    receptor.write_text(donor(1, 22.9) + water_partner)
    waters.write_text("2.9 0 0\n")
    terms, credits, _ = check_report(fixed(binary, receptor, ligand, flags, mode="satisfied_water").stdout)
    require(min(credits) > 0 and same((sum(credits),), (terms[3],)), "Credits were combined across sites")


def check_real(binary, baseline, directory):
    receptor = ROOT / "examples/basic/1iep_receptor.pdbqt"
    ligand = ROOT / "examples/basic/1iep_ligand.pdbqt"
    waters = directory / "waters.xyz"
    # Deliberately distant oxygen for a stable real receptor CLI parity probe.
    waters.write_text("1000 1000 1000\n")
    flags = ("--net_polar_waters", str(waters))
    old = fixed(baseline, receptor, ligand, flags, mode="water").stdout
    new = fixed(binary, receptor, ligand, flags, mode="satisfied_water").stdout
    require(diagnostics(old) == diagnostics(new), "Native energy or old diagnostics changed")
    require(check_report(new)[1][0] > 0, "Real complex has no detected direct credit")
    budget = ("--pose_selection", "late", "--refine_candidates", "6")
    previous_run, previous = dock(baseline, directory, "old-water", budget + ("--net_polar", "water") + flags)
    current_run, current = dock(binary, directory, "new-water", budget + ("--net_polar", "water") + flags)
    previous_run.check_returncode()
    current_run.check_returncode()
    require(previous.read_bytes() == current.read_bytes(), "Legacy global water output changed")
    expected = {}
    for block in filter(str.strip, previous.read_text().split("ENDMDL")):
        native_score, native_rank = values(block, "REMARK VINA NATIVE:")
        expected[native_rank] = (native_score, [line for line in block.splitlines() if line.startswith("ATOM")])
    for index, mode in enumerate(("satisfied", "satisfied_water")):
        options = budget + ("--net_polar", mode) + (flags if index else ())
        run, output = dock(binary, directory, mode, options)
        run.check_returncode()
        require(COUNTS.findall(run.stdout) == COUNTS.findall(previous_run.stdout), "Search budgets changed")
        order = []
        for block in filter(str.strip, output.read_text().split("ENDMDL")):
            native_score, native_rank = values(block, "REMARK VINA NATIVE:")
            atoms = [line for line in block.splitlines() if line.startswith("ATOM")]
            require((native_score, atoms) == expected[native_rank], "Saved native pose or score changed")
            delta = values(block, SATISFACTION_REMARK)[4 + index]
            selected = values(block, "REMARK VINA NET POLAR RESULT:")[0]
            require(same((selected,), (native_score + delta,)), "Wrong selected global score")
            order.append((selected, native_rank))
        require(order == sorted(order) and sorted(rank for _, rank in order) == sorted(expected),
                "Global arm did not sort the complete shared native pool")

    first_model = output.read_text().splitlines()
    require(first_model[0].startswith("MODEL "), "Expected a saved pose")
    remark = "REMARK VINA provenance: retain this input note"
    annotated = directory / "from-saved-satisfied-pose.pdbqt"
    annotated.write_text("\n".join([remark, *first_model[1:first_model.index("ENDMDL")]]) + "\n")
    for mode in ("native", "satisfied"):
        options = budget + (("--net_polar", mode) if mode == "satisfied" else ())
        run, rerun = dock(binary, directory, f"annotated-{mode}", options, ligand=annotated)
        run.check_returncode()
        for block in filter(str.strip, rerun.read_text().split("ENDMDL")):
            require(block.count(SATISFACTION_REMARK) == (mode == "satisfied"),
                    "A saved satisfaction remark leaked into the new pose")
            require(block.count(remark) == 1, "A non-score input remark changed")

    def score_pool(executable, name, options=()):
        output = directory / f"{name}.json"
        subprocess.run([sys.executable, str(ROOT / "tools/score_physical_pool.py"),
                        "--binary", str(executable), "--receptor", str(receptor),
                        "--pool", str(previous), "--output", str(output), *options],
                       capture_output=True, text=True, timeout=30, check=True)
        result = json.loads(output.read_text())
        if result["strain"]:
            del result["strain"]["wall_seconds"]
        for row in result["poses"]:
            del row["engine"]["wall_seconds"]
        return result

    legacy = score_pool(baseline, "legacy-json")
    require(score_pool(binary, "default-json") == legacy, "Default fixed-pool JSON changed")
    dry = score_pool(binary, "dry-json", ("--satisfaction",))
    wet = score_pool(binary, "wet-json", ("--satisfaction", "--waters", str(waters)))
    require(dry["orders"]["satisfied"] == wet["orders"]["satisfied"] == wet["orders"]["satisfied_water"],
            "A far water altered dry/wet ranking")
    require(not dry["orders"]["satisfied_water"] and "satisfied_water" in dry["abstentions"],
            "Fixed-pool water arm failed to abstain")
    for row in dry["poses"]:
        require(row["satisfied_water"] is None and
                same((row["satisfied"],), (row["native"] + row["engine"]["satisfied_delta"],)),
                "Fixed pool did not add the disclosed delta to its saved baseline")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path)
    args = parser.parse_args()
    with TemporaryDirectory(prefix="vina-satisfaction-check-") as directory:
        binary, baseline = args.binary.resolve(strict=True), args.baseline.resolve(strict=True)
        check_analytic(binary, Path(directory))
        check_real(binary, baseline, Path(directory))
    print("Net polar satisfaction: site geometry and credit competition, real 1IEP pool, legacy CLI and JSON passed.")


if __name__ == "__main__":
    main()
