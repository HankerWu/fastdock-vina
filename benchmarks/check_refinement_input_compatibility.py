#!/usr/bin/env python3
"""Exercise pseudoatom geometry and attached-score filtering through Vina and vina_split."""

import argparse
from itertools import combinations
from math import dist
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from check_late_pose_selection import COUNTS, require
from check_pro_rerank import BOX, LIGAND, atoms, models


LATE = ("--pose_selection", "late")
ARMIJO = ("--refinement_optimizer", "armijo")
GAUSS2 = ("--rerank_gauss2_scale", "0.5")
NATIVE_MARKS = ("REMARK VINA RESULT:", "REMARK INTER + INTRA:", "REMARK INTER:",
                "REMARK INTRA:", "REMARK UNBOUND:")
PRO_MARKS = ("REMARK VINA PRO RESULT:", "REMARK VINA NATIVE:", "REMARK VINA PRO GAUSS2 SCALE:",
             *NATIVE_MARKS[1:])
SCORE_MARKS = (*NATIVE_MARKS, *PRO_MARKS[:3], "REMARK CONF_INDEPENDENT:")
USER_REMARKS = ("REMARK VINA provenance: test ligand", "REMARK INTERACTION NOTE: preserve exactly")


def command(args):
    result = subprocess.run(args, check=False, capture_output=True, text=True, timeout=90)
    require(result.returncode == 0, result.stderr or result.stdout)
    return result.stdout


def dock(binary, directory, label, options=(), ligand=LIGAND, seed=1, modes=3, threshold=1):
    output = directory / f"{label}.pdbqt"
    log = command([str(binary), "--ligand", str(ligand), "--receptor", str(BOX.with_suffix("").with_suffix(".pdbqt")),
                   "--config", str(BOX), "--out", str(output), "--cpu", "1", "--seed", str(seed),
                   "--exhaustiveness", "3", "--max_evals", "2000", "--num_modes", str(modes),
                   "--min_rmsd", str(threshold), "--energy_range", "100", *options])
    return output, log


def dummy_check(new, old, directory):
    # A real C and legal W in one rigid root. W may swing farther than C across orientations.
    carbon = "ATOM      1  C   UNL     1      15.190  53.903  16.917  1.00  0.00     0.000 C "
    water_dummy = "ATOM      2  W   UNL     1      17.990  53.903  16.917  1.00  0.00     0.000 W"
    ligand = directory / "carbon-water-dummy.pdbqt"
    ligand.write_text(f"ROOT\n{carbon}\n{water_dummy}\nENDROOT\nTORSDOF 0\n")
    old_native, _ = dock(old, directory, "dummy-old-default", ligand=ligand, modes=9)
    new_native, _ = dock(new, directory, "dummy-new-default", ligand=ligand, modes=9)
    require(old_native.read_bytes() == new_native.read_bytes(), "Default dummy geometry changed")

    def min_real_separation(output):
        coords = [tuple(float(atoms(m)[0][i:i+8]) for i in (30, 38, 46)) for m in models(output)]
        require(len(coords) >= 2, "Dummy fixture did not produce multiple output poses")
        return min(dist(a, b) for a, b in combinations(coords, 2))

    for threshold in (1, 2.5):
        before, before_log = dock(old, directory, f"dummy-old-{threshold}", LATE, ligand=ligand,
                                  modes=9, threshold=threshold)
        after, after_log = dock(new, directory, f"dummy-new-{threshold}", LATE, ligand=ligand,
                                modes=9, threshold=threshold)
        before_min = min_real_separation(before)
        after_min = min_real_separation(after)
        require(before_min < threshold - 0.002 and after_min >= threshold - 0.002,
                f"Dummy fixture did not reproduce and remove subthreshold chemical poses at {threshold} Å")
        (old_counts,) = [tuple(map(int, m)) for m in COUNTS.findall(before_log)]
        (new_counts,) = [tuple(map(int, m)) for m in COUNTS.findall(after_log)]
        require(old_counts[0] == new_counts[0], "Late changed minima collection within trajectories")
        if threshold == 2.5:
            require(old_counts[1] > new_counts[1] > new_counts[3],
                    f"Fixture must exercise both pre-refinement and post-refinement filtering: {old_counts}, {new_counts}")
        print(f"Dummy at {threshold} Å: min C separation {before_min:.3f} → {after_min:.3f} Å, "
              f"counters {old_counts} → {new_counts}")


def split_check(new, old, splitter, directory):
    ligand = directory / "annotated.pdbqt"
    ligand.write_text("\n".join(USER_REMARKS) + "\n" + LIGAND.read_text())
    for origin, origin_flags in (("native", ()), ("pro", GAUSS2)):
        original, _ = dock(old, directory, f"orig-{origin}", origin_flags, ligand=ligand, modes=1)
        prefix = directory / f"split-{origin}-"
        command([str(splitter), "--input", str(original), "--ligand", str(prefix)])
        (split,) = directory.glob(f"{prefix.name}*.pdbqt")
        original_marks = PRO_MARKS if origin == "pro" else NATIVE_MARKS
        require(tuple(line.partition(":")[0] + ":" for line in split.read_text().splitlines()
                      if line.startswith(SCORE_MARKS)) == original_marks, "Split input lacks its old scores")

        legacy, _ = dock(old, directory, f"default-old-{origin}", ligand=split, modes=1)
        default, _ = dock(new, directory, f"default-new-{origin}", ligand=split, modes=1)
        require(legacy.read_bytes() == default.read_bytes(), "Default output changed for split input")

        for name, flags in (("late", LATE), ("armijo", ARMIJO)):
            stale, _ = dock(old, directory, f"stale-{origin}-{name}", flags, ligand=split, modes=1)
            (stale_model,) = models(stale)
            require(sum(line.startswith(SCORE_MARKS) for line in stale_model) > len(NATIVE_MARKS),
                    f"Optimizer baseline did not reproduce stale marks for {origin}/{name}")
            clean, _ = dock(new, directory, f"clean-{origin}-{name}", flags, ligand=split, modes=1)
            (clean_model,) = models(clean)
            require(tuple(line.partition(":")[0] + ":" for line in clean_model
                          if line.startswith(SCORE_MARKS)) == NATIVE_MARKS,
                    f"Wrong score headers for {origin}/{name}")
            require(all(clean_model.count(line) == 1 for line in USER_REMARKS),
                    f"User remarks changed for {origin}/{name}")
            require([line for line in clean_model if not line.startswith(SCORE_MARKS)] ==
                    [line for line in stale_model if not line.startswith(SCORE_MARKS)],
                    f"Non-score lines changed for {origin}/{name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("splitter", type=Path)
    args = parser.parse_args()
    new, old, splitter = (arg.resolve(strict=True) for arg in (args.binary, args.baseline, args.splitter))
    with TemporaryDirectory(prefix="fastdock-refinement-smoke-") as work:
        directory = Path(work)
        for label, flags in (("native", ()), ("late", LATE), ("armijo", ARMIJO),
                             ("all", LATE + ARMIJO + GAUSS2)):
            before, _ = dock(old, directory, f"regular-old-{label}", flags, seed=42)
            after, _ = dock(new, directory, f"regular-new-{label}", flags, seed=42)
            require(before.read_bytes() == after.read_bytes(),
                    f"Clean non-dummy output differs from optimizer baseline for {label}")
        dummy_check(new, old, directory)
        split_check(new, old, splitter, directory)
    print("Clean-input byte parity, both late geometry filters, and real split score-header regressions passed.")


if __name__ == "__main__":
    main()
