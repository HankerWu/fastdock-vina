#!/usr/bin/env python3
"""Check trajectory archive CLI output, ranking and limits on the 1IEP fixture."""

import argparse
import math
from pathlib import Path
import re
import subprocess
import tempfile
from typing import NamedTuple

from check_late_pose_selection import COUNTS, ROOT, require


class Pose(NamedTuple):
    atoms: tuple[str, ...]
    score: float
    rank: int
    source: str
    votes: int
    eligible: bool


HEADERS = ("ARCHIVE RESULT", "NATIVE", "ARCHIVE SUPPORT", "ARCHIVE SOURCE", "ARCHIVE MODE")
STATS = re.compile(r"Trajectory archive: mode=(\w+) ([\w= ]+)")


def parse_output(text, mode):
    poses = []
    for block in text.split("ENDMDL"):
        if not block.strip():
            continue
        lines = block.strip().splitlines()
        fields = []
        for header, line in zip(HEADERS, lines[1:6]):
            prefix = f"REMARK VINA {header}: "
            require(line.startswith(prefix) and sum(x.startswith(prefix) for x in lines) == 1,
                    f"Missing or repeated {header}")
            fields.append(line.removeprefix(prefix).split())
        result, native, support, source, archive_mode = fields
        votes, eligible, computed = map(int, support)
        require(result[0] == native[0] and all(math.isfinite(float(x)) for x in result)
                and archive_mode == [mode] and 0 <= votes <= 6 and eligible in (0, 1)
                and computed == (mode == "support") and (computed or votes == 0)
                and source in (["historical"], ["accepted"]), "Invalid archive metadata")
        atoms = tuple(line for line in lines if line.startswith(("ATOM  ", "HETATM")))
        require(atoms, "Empty ligand")
        poses.append(Pose(atoms, float(native[0]), int(native[1]), source[0], votes, bool(eligible)))
    require(poses, "No poses")
    return poses


def check_ranking(accepted, support):
    native_poses = parse_output(accepted, "accepted")
    support_poses = parse_output(support, "support")
    require(native_poses == sorted(native_poses, key=lambda p: p.score)
            and [p.rank for p in native_poses] == list(range(1, len(native_poses) + 1)),
            "Native scores or ranks are out of order")
    bound = native_poses[0].score + 1
    require(all(p.eligible == (p.score <= bound) for p in native_poses + support_poses),
            "Eligibility differs from the native +1 kcal/mol window")
    require(support_poses == sorted(support_poses,
                                   key=lambda p: (not p.eligible, -p.votes if p.eligible else 0, p.rank)),
            "Support output contradicts its recorded ranking")
    require(native_poses == [p._replace(votes=0) for p in sorted(support_poses, key=lambda p: p.rank)],
            "The modes disagree on native scores, ranks, coordinates, eligibility or source")
    return native_poses, support_poses


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    baseline = args.baseline.resolve(strict=True)
    fixture = ROOT / "examples/basic"
    ligand = fixture / "1iep_ligand.pdbqt"
    late = ("--pose_selection", "late")

    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)

        def run(name, options=(), engine=binary, cpu=1, range_kcal=1e9, ligand_file=ligand, seed=3, error=None):
            output = directory / f"{name}.pdbqt"
            cmd = [str(engine), "--receptor", str(fixture / "1iep_receptor.pdbqt"),
                   "--ligand", str(ligand_file), "--config", str(fixture / "1iep_receptor.box.txt"),
                   "--out", str(output), "--cpu", str(cpu), "--seed", str(seed), "--exhaustiveness", "6",
                   "--max_evals", "1100", "--num_modes", "6", "--energy_range", str(range_kcal)]
            result = subprocess.run(cmd + list(options), text=True, capture_output=True, timeout=90)
            if error:
                require(result.returncode != 0 and error in result.stderr and not output.exists(), result.stderr)
                return None
            require(result.returncode == 0, result.stderr)
            return result.stdout, output.read_text()

        for name, options in (("default", ()), ("late", late)):
            _, previous = run("previous-" + name, options, engine=baseline)
            log, current = run("current-" + name, options)
            require(previous == current and not STATS.search(log), "Disabled archive changed byte output")

        def archive_run(name, mode, options=(), **kwargs):
            log, contents = run(name, late + ("--trajectory_archive", mode) + options, **kwargs)
            ((actual_mode, stats),) = STATS.findall(log)
            require(actual_mode == mode and "Late pose selection:" in log, "Unexpected archive accounting")
            return dict((key, int(value)) for key, value in (item.split("=") for item in stats.split())), contents

        stats, accepted = archive_run("accepted", "accepted")
        support_stats, support = archive_run("support", "support")
        require(stats == support_stats, "Ranking changed the common refinement pool")
        require(stats["trajectories"] == 6 and 0 < stats["native"] <= 36 and 0 < stats["extra"] <= 36
                and 0 < stats["extra_prescore"] <= stats["extra"]
                and stats["prescore"] == stats["native"] + stats["extra_prescore"]
                and 0 < stats["extra_refine"] <= stats["refine"] <= 6
                and stats["refine"] <= stats["refinement_calls"] <= 5 * stats["refine"]
                and stats["native_search_evals"] >= 6 * 1100 and stats["native_search_refinements"] > 0,
                f"Incorrect archive accounting or candidate bounds {stats}")
        native_poses, support_poses = check_ranking(accepted, support)
        require(len(native_poses) == len(support_poses) == stats["retained"] <= stats["refine"],
                "Wrong refined-pose cap")
        require(any(p.source == "accepted" for p in native_poses), "The fixture did not retain an accepted seed")

        threaded_stats, threaded = archive_run("threaded", "support", cpu=3)
        require(threaded_stats == stats and threaded == support, "Support output depends on worker count")
        narrow_stats, narrow = archive_run("narrow", "support", range_kcal=0.001)
        require(narrow_stats == stats and parse_output(narrow, "support") ==
                [p for p in support_poses if p.score <= native_poses[0].score + 0.001],
                "Reported poses do not follow the native energy window")

        same_cap = ("--refine_candidates", "6")
        for mode, expected in (("accepted", accepted), ("support", support)):
            same_stats, same = archive_run("same-" + mode, mode, same_cap)
            require(same_stats == stats and same == expected, "Explicit default cap changed archive output")
        late_default_log, late_default = run("late-default", late)
        _, late_same = run("late-same", late + same_cap)
        require(late_default == late_same, "Explicit default cap changed native Late output")
        large_cap = ("--refine_candidates", "18")
        late_large_log, _ = run("late-large", late + large_cap)
        (late_default_counts,) = COUNTS.findall(late_default_log)
        (late_large_counts,) = COUNTS.findall(late_large_log)
        saved, refined, calls, retained = map(int, late_large_counts)
        require(saved == int(late_default_counts[0]) and 6 < refined <= 18 and refined <= calls <= 5 * refined
                and 0 < retained <= 6, "Expanded native Late changed its search pool or output limit")
        old_stats, _ = archive_run("default-budget", "support", seed=5)
        large_stats, expanded = archive_run("large-support", "support", large_cap, seed=5)
        accepted_stats, expanded_accepted = archive_run("large-accepted", "accepted", large_cap, seed=5)
        for key in ("trajectories", "native", "extra", "prescore", "extra_prescore",
                    "native_search_evals", "native_search_refinements"):
            require(old_stats[key] == large_stats[key], "Higher refinement cap altered upstream search")
        require(large_stats == accepted_stats and old_stats["refine"] == 6 and large_stats["refine"] == 18
                and old_stats["retained"] < large_stats["retained"] <= 6
                and large_stats["refine"] <= large_stats["refinement_calls"] <= 5 * large_stats["refine"],
                f"Expanded cap did not recover distinct poses with bounded retries {old_stats}, {large_stats}")
        expanded_native, expanded_support = check_ranking(expanded_accepted, expanded)
        require(len(expanded_support) == large_stats["retained"], "Expanded pool exceeds its output cap")
        _, narrow_expanded = archive_run("narrow-expanded", "support", large_cap, seed=5, range_kcal=0.001)
        require(parse_output(narrow_expanded, "support") ==
                [p for p in expanded_support if p.score <= expanded_native[0].score + 0.001],
                "Expanded archive did not honor the native energy window")

        previous_pose = directory / "previous-pose.pdbqt"
        previous_pose.write_text(support.split("ENDMDL", 1)[0].split("\n", 1)[1])
        _, redocked = archive_run("redocked", "accepted", ligand_file=previous_pose)
        parse_output(redocked, "accepted")

        archive = late + ("--trajectory_archive", "accepted")
        for label, options, message in (
            ("unknown", late + ("--trajectory_archive", "unknown"), "trajectory_archive must be accepted or support"),
            ("native", ("--trajectory_archive", "accepted"), "trajectory_archive requires"),
            ("armijo", archive + ("--refinement_optimizer", "armijo"), "trajectory_archive requires"),
            ("polar", archive + ("--polar_geometry", "distance"), "trajectory_archive requires"),
            ("template", archive + ("--pose_template", "unused", "--template_stage", "rank"),
             "trajectory_archive requires"),
            ("gaussian", archive + ("--rerank_gauss2_scale", "0.5"), "trajectory_archive requires"),
            ("local", archive + ("--local_only",), "late pose selection requires"),
            ("multiple", archive + ("--ligand", str(ligand)), "late pose selection requires"),
            ("small-cap", late + ("--refine_candidates", "5"), "refine_candidates must be positive and at least"),
            ("zero-cap", late + ("--refine_candidates", "0"), "refine_candidates must be positive and at least"),
            ("cap-without-late", large_cap, "refine_candidates requires Late Vina"),
            ("cap-with-armijo", late + large_cap + ("--refinement_optimizer", "armijo"), "refine_candidates requires"),
            ("cap-with-polar", late + large_cap + ("--polar_geometry", "distance"), "refine_candidates requires"),
        ):
            run(label, options, error=message)

        for atom_type in ("G0", "W", "CG0"):
            lines = ligand.read_text().splitlines()
            index = next(i for i, line in enumerate(lines) if line.startswith("ATOM  "))
            lines[index] = lines[index][:77] + atom_type
            unsupported = directory / f"{atom_type}.pdbqt"
            unsupported.write_text("\n".join(lines) + "\n")
            run("bad-" + atom_type, archive, ligand_file=unsupported,
                error="trajectory_archive does not support ligand dummy atoms or macrocycle closure types")

    print("Disabled parity, candidate caps, provenance, ranking, energy window and option validation passed.")


if __name__ == "__main__":
    main()
