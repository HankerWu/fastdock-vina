"""Paired, fixed-budget X-ray redocking and cross-docking of a graph-distance fix."""

import argparse
import csv
import hashlib
import json
import random
import sys
from pathlib import Path

CASES = {
    "3DCT": ("NR1H4", "064", "A"),
    "3FLI": ("NR1H4", "33Y", "A"),
    "2XP2": ("ALK", "VGH", "A"),
    "3K54": ("BTK", "1N1", "A"),
    "1CBS": ("CRABP2", "REA", "A"),
    "6WTO": ("JAK2", "3JW", "A"),
    "6VN8": ("JAK2", "3JW", "B"),
    "4W9X": ("BMP2K", "3JW", "A"),
    "6VGL": ("JAK2", "RXT", "D"),
    "5I3O": ("BMP2K", "IDV", "A"),
    "6IG8": ("CSF1R", "A7O", "A"),
}
SEED = 817139
CAPS = (3000, 12000)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def total_evals(row):
    return sum(int(row[s + "_evals"]) for s in ("hunt", "second", "explicit"))


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument(
        "action", choices=["prepare", "run", "cost", "replay", "confirm", "production"]
    )
    a.add_argument("workspace", type=Path)
    a.add_argument("--repo", required=True, type=Path)
    a.add_argument("--upstream", type=Path)
    a.add_argument("--counts", type=Path)
    a.add_argument("--parent", type=Path)
    a.add_argument("--fixed", type=Path)
    a.add_argument("--release", type=Path)
    a.add_argument("--ids", nargs="*")
    args = a.parse_args()
    work = args.workspace.resolve()
    sys.path.insert(0, str(args.repo.resolve() / "benchmarks/search_backtracking/pose"))
    import backtracking_panel as p

    p.CASES = CASES
    # The generic receptor preparer retains deposited magnesium (e.g. CSF1R).
    p.paired.RADII["Mg"] = 1.73
    for i in range(4):
        p.paired.RADII["CG" + str(i)] = p.paired.RADII["C"]
    original_atoms = p.paired.atoms

    def real_atoms(lines):
        return original_atoms(
            [
                l
                for l in lines
                if not l.startswith(("ATOM  ", "HETATM"))
                or l.split()[-1] not in {"G0", "G1", "G2", "G3"}
            ]
        )

    p.paired.atoms = real_atoms
    if args.action == "prepare":
        for pdb in args.ids or CASES:
            p.prepare_input(work, pdb)
        return
    data = {q: json.loads((work / q / "preparation.json").read_text()) for q in CASES}
    # Map the depositors' two numbering conventions into the same UniProt
    # Q96RI1 coordinate system using each PDB's deposited DBREF record.
    for q in ("3FLI", "3DCT"):
        row = next(
            x.split()
            for x in (work / "raw" / (q + ".pdb")).read_text().splitlines()
            if x.startswith("DBREF") and "Q96RI1" in x
        )
        shift = int(row[-2]) - int(row[3])
        if shift != int(row[-1]) - int(row[4]):
            raise ValueError("Nonlinear NR1H4 residue numbering")
        data[q]["ca"] = {
            str(int(k.split(":")[0]) + shift) + ":" + k.split(":")[1]: v
            for k, v in data[q]["ca"].items()
        }
    if args.action == "production":
        rows = []
        for summary in ("cells.csv", "confirmation.csv"):
            with (work / summary).open() as handle:
                baseline = [r for r in csv.DictReader(handle) if r["arm"] == "fixed"]
            for old in baseline:
                row, _ = p.paired.execute(
                    work,
                    data,
                    old["query"],
                    old["receiver"],
                    int(old["seed"]),
                    int(old["cap"]),
                    "production",
                    args.fixed,
                    False,
                )
                new_digest = digest(work / "runs" / (row["key"] + ".pdbqt"))
                if old["sha256"] != new_digest or any(
                    str(row[s + "_" + k]) != old[s + "_" + k]
                    for s in ("hunt", "second", "explicit")
                    for k in ("calls", "evals")
                ):
                    raise ValueError(
                        "Integrated build differs from isolated BFS: " + row["key"]
                    )
                rows.append(
                    {
                        "key": row["key"],
                        "sha256": new_digest,
                        "seconds": row["seconds"],
                        "binary_sha256": digest(args.fixed),
                    }
                )
        p.paired.write(work / "production.csv", rows)
        return
    if args.action == "confirm":
        rows, modes = [], []
        with (work / "cells.csv").open() as handle:
            eligible = [
                r
                for r in csv.DictReader(handle)
                if r["arm"] == "parent" and CASES[r["query"]][1] == "3JW"
            ]
        random.Random(191009).shuffle(eligible)
        for i, previous in enumerate(eligible):
            for arm in ("fixed", "parent") if i % 2 else ("parent", "fixed"):
                row, poses = p.paired.execute(
                    work,
                    data,
                    previous["query"],
                    previous["receiver"],
                    191009,
                    int(previous["cap"]),
                    arm,
                    getattr(args, arm),
                    False,
                )
                row["sha256"] = digest(work / "runs" / (row["key"] + ".pdbqt"))
                rows.append(row)
                modes.extend(poses)
        p.paired.write(work / "confirmation.csv", rows)
        p.paired.write(work / "confirmation_modes.csv", modes)
        return
    if args.action in ("cost", "replay"):
        with (work / "cells.csv").open() as f:
            main_rows = list(csv.DictReader(f))
        index = {(r["query"], r["receiver"], r["cap"], r["arm"]): r for r in main_rows}
        rows, modes = [], []
        for old in main_rows:
            q, r, c, arm = old["query"], old["receiver"], old["cap"], old["arm"]
            if CASES[q][1] != "3JW":
                continue
            if args.action == "cost":
                if arm != "fixed" or total_evals(old) <= total_evals(
                    index[q, r, c, "parent"]
                ):
                    continue
                requested = int(int(c) * 0.85)
                row, poses = p.paired.execute(
                    work, data, q, r, SEED, requested, "cost", args.fixed, False
                )
                row["parent_cap"] = c
                rows.append(row)
                modes.extend(poses)
            elif arm == "upstream":
                row, _ = p.paired.execute(
                    work, data, q, r, SEED, int(c), "release", args.release, False
                )
                identical = (
                    digest(work / "runs" / (row["key"] + ".pdbqt")) == old["sha256"]
                )
                rows.append(
                    {
                        "query": q,
                        "receiver": r,
                        "cap": c,
                        "identical": identical,
                        "seconds": row["seconds"],
                        "binary_sha256": digest(args.release),
                    }
                )
                if not identical:
                    raise ValueError("Official release differed at " + row["key"])
        if rows:
            p.paired.write(work / (args.action + ".csv"), rows)
        if modes:
            p.paired.write(work / "cost_modes.csv", modes)
        return
    p.write_sources(work)
    (work / "runs").mkdir(exist_ok=True)
    cells = [
        (q, r, c)
        for q in CASES
        for r in CASES
        if CASES[q][0] == CASES[r][0]
        and data[r]["receptor_status"] == "ready"
        and (q == r or CASES[q][1] != CASES[r][1] and q != "6WTO")
        for c in CAPS
    ]
    eligible, rejected = [], []
    for q, r, c in cells:
        try:
            p.paired.base.alignment(data[q], data[r])
        except ValueError as err:
            rejected.append({"query": q, "receiver": r, "cap": c, "reason": str(err)})
        else:
            eligible.append((q, r, c))
    cells = eligible
    if rejected:
        p.paired.write(work / "alignment_exclusions.csv", rejected)
    random.Random(SEED).shuffle(cells)
    rows, modes = [], []
    for i, (q, r, c) in enumerate(cells):
        arms = ["upstream", "parent", "fixed"]
        for arm in arms[i % 3 :] + arms[: i % 3]:
            row, poses = p.paired.execute(
                work, data, q, r, SEED, c, arm, getattr(args, arm), False, args.counts
            )
            row["sha256"] = digest(work / "runs" / (row["key"] + ".pdbqt"))
            if any(
                row[k + "_evals"] == "" or row[k + "_calls"] == ""
                for k in ("hunt", "second", "explicit")
            ):
                raise ValueError("Missing counter for " + row["key"])
            rows.append(row)
            modes.extend(poses)
    p.paired.write(work / "cells.csv", rows)
    p.paired.write(work / "modes.csv", modes)
    (work / "binaries.json").write_text(
        json.dumps(
            {
                k: digest(getattr(args, k))
                for k in ("upstream", "counts", "parent", "fixed")
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
