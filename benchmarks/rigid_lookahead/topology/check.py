"""Check every blind pre-refinement start, preserving results separately by panel/ligand/arm."""

import argparse
import csv
import itertools
import subprocess
from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory


def positive_control(probe, ligand):
    # Pull only the deposited-input model's explicit donor hydrogen away from N.
    # This corrupts an artificial copy and must change donor typing and the graph.
    with TemporaryDirectory() as scratch:
        changed = Path(scratch) / "detached_H.pdbqt"
        lines = ligand.read_text().splitlines()
        h = next(i for i, line in enumerate(lines) if line.startswith("ATOM") and
                 line[6:11].strip() == "13" and line.split()[-1] == "HD")
        lines[h] = lines[h][:30] + f"{90.0:8.3f}" + lines[h][38:]
        changed.write_text("\n".join(lines) + "\n")
        out = subprocess.check_output([str(probe)], text=True, input=f"{ligand}\t{changed}\n")
        if out.split()[:2] != ["1", "1"]:
            raise ValueError("Failed to detect the deliberate donor/edge change")


def main(a):
    positive_control(a.probe, a.dev / "5EPS" / "ligand.pdbqt")
    sets = [(f"{stage}_{budget}_{arm}", getattr(a, stage), a.sweep / f"{stage}_{budget}_{arm}")
            for stage, budget, arm in itertools.product(("dev", "secondary"), (9, 24, 48), ("raw", "sham", "gradient"))]
    sets.append(("RBP4_9", a.holdout, a.rbp4_runs))
    jobs = []
    for label, base, run in sets:
        # All starts enter refinement regardless of subsequent ranking / deduplication.
        starts = sorted((run / "runs").glob("*/*.pdbqt"))
        with (run / "poses.csv").open() as f:
            rows = list(csv.DictReader(f))
        expected = {run / "runs" / "_".join(row[k] for k in ("receiver", "query", "seed", "arm"))
                    / (row["slot"] + ".pdbqt") for row in rows}
        if expected != set(starts) or len(starts) != len(rows) or not starts:
            raise ValueError("Missing, extra, or duplicated blind starts: " + label)
        for path in starts:
            _, query, _, arm = path.parent.name.split("_")
            jobs.append((label, query, arm, base / query / "ligand.pdbqt", path))
    output = subprocess.check_output([str(a.probe)], text=True, input="".join(
        f"{original}\t{path}\n" for _, _, _, original, path in jobs))
    lines = output.splitlines()
    if len(lines) != len(jobs):
        raise ValueError("Incomplete parser response")
    stats = defaultdict(lambda: [0] * 7)
    mismatches = []
    for (label, query, arm, _, path), line in zip(jobs, lines):
        *values, delta = line.split()
        flags = list(map(int, values))
        if len(flags) != 5:
            raise ValueError("Wrong parser result: " + line)
        row = stats[label, query, arm]
        row[0] += 1
        for i, flag in enumerate(flags, 1):
            row[i] += flag
        row[-1] = max(row[-1], float(delta))
        if any(flags):
            mismatches.append([label, query, arm, path.parent.name, path.name, *flags])
    names = ["atom_type_or_charge", "covalent_edges", "rotatable_edges", "branch_or_dof", "pair_or_table_index"]
    a.output.mkdir(parents=True, exist_ok=True)
    with (a.output / "summary.csv").open("w") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["panel", "query", "arm", "starts", *names, "largest_bond_length_change_A"])
        for group, row in sorted(stats.items()):
            w.writerow([*group, *row])
    with (a.output / "mismatches.csv").open("w") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["panel", "query", "arm", "cell", "file", *names])
        w.writerows(mismatches)
    print(f"{len(jobs)} starts, {len(mismatches)} with a changed topology")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("probe", "dev", "secondary", "sweep", "holdout", "rbp4_runs", "output"):
        p.add_argument(arg, type=Path)
    main(p.parse_args())
