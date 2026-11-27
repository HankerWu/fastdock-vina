"""Reproduce fixed-budget native-first-pose comparisons from frozen input fixtures."""

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
import re
import subprocess
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
CAP = 12000
SEED = 614039
CAL_SEED = 614033
COUNT = (
    "hunt_evals",
    "second_evals",
    "explicit_evals",
    "hunt_calls",
    "second_calls",
    "explicit_calls",
)
RECS = {"Gal3": ("1KJL", "3T1L"), "CAII": ("6WQ4", "6WQ8")}
CA_LIGS = ("6WQ4", "6WQ5", "6WQ7", "6WQ8", "6WQ9")


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w", newline="") as stream:
        out = csv.DictWriter(stream, list(rows[0]), lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def cells(work, panel):
    if panel == "Gal3":
        fixtures = work / "gal_inputs"
        if not fixtures.exists():
            with tarfile.open(
                BENCH / "pair_gradient/galectin_external/raw.tar.xz"
            ) as ar:
                ar.extractall(
                    work / "gal",
                    members=[m for m in ar if m.name.startswith("prepared/")],
                    filter="data",
                )
            (work / "gal/prepared").rename(fixtures)
        compounds = sorted(
            (
                r["compound"]
                for r in read(BENCH / "pair_gradient/galectin_external/manifest.csv")
            ),
            key=int,
        )
        if len(set(compounds)) != 47:
            raise ValueError("Incomplete Gal3 panel")
        return [
            (
                panel,
                rec,
                c,
                fixtures / (rec + ".pdbqt"),
                fixtures / (rec + ".box"),
                fixtures / (c + ".pdbqt"),
            )
            for rec, c in itertools.product(RECS[panel], compounds)
        ]
    fixtures = work / "ca/prepared"
    return [
        (
            panel,
            rec,
            c,
            fixtures / (rec + "_receptor.pdbqt"),
            fixtures / (rec + "_box.txt"),
            fixtures / (c + "_ligand.pdbqt"),
        )
        for rec, c in itertools.product(RECS[panel], CA_LIGS)
    ]


def dock(cell, work, binary, arm, seed, cap, include_counters, diag=False):
    panel, rec, compound, protein, box, ligand = cell
    key = f"{panel}_{rec}_{compound}_{seed}_{arm}_{cap}"
    raw = work / "raw"
    outfile = raw / (key + (".replay" if diag else ".pdbqt"))
    args = [
        binary,
        "--receptor",
        protein,
        "--config",
        box,
        "--ligand",
        ligand,
        "--seed",
        str(seed),
        "--cpu",
        "1",
        "--exhaustiveness",
        "1",
        "--max_evals",
        str(cap),
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--verbosity",
        "0",
        "--out",
        outfile,
    ]
    if include_counters:
        args.append("--search_counters")
    start = time.perf_counter_ns()
    p = subprocess.run(args, check=True, text=True, capture_output=True)
    elapsed = time.perf_counter_ns() - start
    log = p.stdout + p.stderr
    (raw / (key + (".counts" if diag else ".log"))).write_text(log)
    all_scores = re.findall(
        r"^REMARK VINA RESULT:\s+([-\d.]+)", outfile.read_text(), re.MULTILINE
    )
    if not all_scores:
        raise ValueError(key)
    values = {k: next(iter(re.findall(r"\b" + k + r"=(\d+)", log)), "") for k in COUNT}
    join = (
        sum(int(values[k]) for k in COUNT[:3])
        if all(values[k] != "" for k in COUNT[:3])
        else ""
    )
    return dict(
        key=key,
        panel=panel,
        receptor=rec,
        compound=compound,
        seed=seed,
        arm=arm,
        cap=cap,
        elapsed_ns=elapsed,
        first_score=all_scores[0],
        modes=len(all_scores),
        joint=join,
        **values,
    )


def run(args):
    work = args.work.resolve()
    (work / "raw").mkdir(parents=True, exist_ok=True)
    os.sched_setaffinity(0, {args.cpu})
    entries = cells(work, args.panel)
    random.Random(SEED).shuffle(entries)
    arms = {"parent": args.parent, "retention": args.candidate}
    name = args.panel + "_" + args.stage
    result = []
    if args.stage == "upstream":
        arms = {"upstream": args.upstream, "upstream_counts": args.upstream_counts}
        for cell in entries:
            row = dock(cell, work, args.upstream, "upstream", SEED, CAP, False)
            counted = dock(
                cell, work, args.upstream_counts, "upstream", SEED, CAP, False, True
            )
            p = work / "raw" / row["key"]
            if (
                p.with_suffix(".pdbqt").read_bytes()
                != p.with_suffix(".replay").read_bytes()
            ):
                raise ValueError("Upstream count replay differs: " + row["key"])
            p.with_suffix(".replay").unlink()
            if not counted["joint"]:
                raise ValueError("Missing upstream counts")
            row.update({k: counted[k] for k in (*COUNT, "joint")})
            result.append(row)
    elif args.stage == "fixed":
        if args.ungated:
            arms["ungated"] = args.ungated
        for i, cell in enumerate(entries):
            order = list(arms)
            j = i % len(order)
            for arm in order[j:] + order[:j]:
                row = dock(cell, work, arms[arm], arm, SEED, CAP, True)
                if not row["joint"]:
                    raise ValueError("Missing actual counters")
                result.append(row)
    elif args.stage == "matched":
        fixed = read(work / (args.panel + "_fixed.csv"))

        def total(arm, col):
            return sum(int(r[col]) for r in fixed if r["arm"] == arm)

        excess = any(
            total("retention", k) > 1.05 * total("parent", k)
            for k in ("joint", "elapsed_ns")
        )
        if not excess:
            print(f"{args.panel}: neither candidate cost exceeds 105% of parent")
            return
        # Intersect the time and joint-call ceilings by linear interpolation.
        calibration = []
        budgets = []
        for i, cell in enumerate(entries):
            specs = [("parent", CAP), ("retention", 3000), ("retention", CAP)]
            j = i % len(specs)
            timing = {}
            for arm, cap in specs[j:] + specs[:j]:
                r = dock(cell, work, arms[arm], arm, CAL_SEED, cap, True)
                calibration.append(r)
                timing[arm, cap] = r
            ceilings = []
            for k in ("elapsed_ns", "joint"):
                t, lo, hi = (
                    int(timing[a, c][k])
                    for a, c in (
                        ("parent", CAP),
                        ("retention", 3000),
                        ("retention", CAP),
                    )
                )
                ceilings.append(
                    CAP
                    if hi <= t
                    else 3000
                    if lo >= t or hi <= lo
                    else 3000 + 9000 * (t - lo) // (hi - lo)
                )
            budgets.append(
                {
                    "panel": cell[0],
                    "receptor": cell[1],
                    "compound": cell[2],
                    "cap": min(ceilings),
                }
            )
        write(work / (args.panel + "_calibration.csv"), calibration)
        write(work / (args.panel + "_budgets.csv"), budgets)
        for cell, b in zip(entries, budgets, strict=True):
            cached = next(
                r
                for r in fixed
                if r["receptor"] == cell[1]
                and r["compound"] == cell[2]
                and r["arm"] == "retention"
            )
            if b["cap"] == CAP:
                result.append(cached)
            else:
                result.append(
                    dock(cell, work, args.candidate, "retention", SEED, b["cap"], True)
                )
    write(work / (name + ".csv"), result)
    (work / (name + "_binary.json")).write_text(
        json.dumps(
            {
                "cpu": sorted(os.sched_getaffinity(0)),
                "binary_sha256": {
                    a: hashlib.sha256(Path(p).read_bytes()).hexdigest()
                    for a, p in arms.items()
                },
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=("upstream", "fixed", "matched"))
    p.add_argument("panel", choices=tuple(RECS))
    p.add_argument("work", type=Path)
    p.add_argument("--cpu", type=int, default=7)
    for opt in ("upstream", "upstream-counts", "parent", "candidate", "ungated"):
        p.add_argument("--" + opt, type=Path)
    run(p.parse_args())
