"""Dock the Sörme Gal3 panel and compare within-assay rank credits on first poses."""

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal
import hashlib
import json
import fcntl
from contextlib import nullcontext
from itertools import combinations
import os
from pathlib import Path
import random
import re
import subprocess
import sys

from scipy.stats import kendalltau, spearmanr

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH / "rigid_flip"))
from holdout import paired

ARMS = (
    "parent",
    "upstream",
    "always",
    "adaptive",
    "always_reduced",
    "adaptive_reduced",
)
RECEPTORS = ("1KJL", "3T1L")


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def settings(args, entries):
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    names = [r["compound"] + ".pdbqt" for r in entries]
    names.extend(r + suffix for r in RECEPTORS for suffix in (".pdbqt", ".box"))
    engines = (
        ["parent", "always", "adaptive"]
        if args.confirm
        else ["parent", "always", "adaptive", "upstream", "counts"]
    )
    return dict(
        cap=args.cap,
        seed=args.seed,
        confirm=args.confirm,
        cpus=args.cpus,
        inputs={name: digest(args.inputs / name) for name in names},
        engines={arm: digest(getattr(args, arm)) for arm in engines},
    )


def lane(args, jobs, cpu):
    os.sched_setaffinity(0, {cpu})
    rows = []
    for receptor, compound in jobs:
        order = ["parent", "always", "adaptive"] if args.confirm else list(ARMS)
        random.Random(71437 + int(compound) + RECEPTORS.index(receptor)).shuffle(order)
        for arm in order:
            name = arm.split("_")[0]
            cap = args.cap * 85 // 100 if "reduced" in arm else args.cap
            stem = args.out / "runs" / f"{receptor}_{compound}_{arm}"
            output = stem.with_suffix(".pdbqt")
            log = stem.with_suffix(".log")
            cost = stem.with_suffix(".cost")
            flags = [
                "--config",
                str(args.inputs / (receptor + ".box")),
                "--receptor",
                str(args.inputs / (receptor + ".pdbqt")),
                "--ligand",
                str(args.inputs / (compound + ".pdbqt")),
                "--cpu",
                "1",
                "--exhaustiveness",
                "1",
                "--seed",
                str(args.seed),
                "--max_evals",
                str(cap),
                "--num_modes",
                "9",
                "--energy_range",
                "100",
            ]
            countlog = stem.with_suffix(".counts")
            # The upstream all-pairs tables can exceed a gigabyte on the bivalent
            # ligands. Serialize upstream processes; acquire the lock before timing.
            lock = (
                (args.out / ".upstream.lock").open("a")
                if name == "upstream"
                else nullcontext()
            )
            with lock as file:
                if name == "upstream":
                    fcntl.flock(file, fcntl.LOCK_EX)
                if not (output.exists() and cost.exists() and log.exists()):
                    subprocess.run(
                        [
                            sys.executable,
                            "-S",
                            str(BENCH / "rigid_flip/process_cost.py"),
                            str(log),
                            str(cost),
                            str(getattr(args, name)),
                            *flags,
                            "--out",
                            str(output),
                            *([] if name == "upstream" else ["--search_counters"]),
                        ],
                        check=True,
                    )
                observed = log.read_text()
                if name == "upstream":
                    if not countlog.exists():
                        replica = stem.with_suffix(".check")
                        result = subprocess.run(
                            [str(args.counts), *flags, "--out", str(replica)],
                            check=True,
                            capture_output=True,
                            text=True,
                        )
                        if replica.read_bytes() != output.read_bytes():
                            raise ValueError(
                                "Upstream counter replay changed poses: " + str(stem)
                            )
                        replica.unlink()
                        countlog.write_text(result.stdout + result.stderr)
                    observed = countlog.read_text()
            stats = paired.counts(observed)
            if any(value == "" for value in stats.values()):
                raise ValueError("Missing counters: " + str(stem))
            torsion = dict(
                re.findall(r"\b(axes|other|offered|redirected)=(\d+)", observed)
            )
            rows.append(
                dict(
                    arm=arm,
                    receptor=receptor,
                    compound=compound,
                    cap=cap,
                    cpu=cpu,
                    score=re.search(
                        r"REMARK VINA RESULT:\s+([-\d.]+)", output.read_text()
                    )[1],
                    sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                    **stats,
                    **json.loads(cost.read_text()),
                    **{
                        k: torsion.get(k, "")
                        for k in ("axes", "other", "offered", "redirected")
                    },
                )
            )
    return rows


def summarize(out, confirm=False):
    rows = read(out / "cells.csv")
    manifest = {
        r["compound"]: r
        for r in read(BENCH / "pair_gradient/galectin_external/manifest.csv")
    }
    arms = ("parent", "always", "adaptive") if confirm else ARMS
    score = {
        (r["arm"], r["receptor"], r["compound"]): Decimal(r["score"]) for r in rows
    }
    assert len(score) == len(arms) * 2 * len(manifest) == len(rows)
    for arm in arms:
        for i in manifest:
            score[arm, "mean", i] = sum(score[arm, rec, i] for rec in RECEPTORS) / 2
    summary, pairs = [], []
    for receptor in (*RECEPTORS, "mean"):
        for arm in arms:
            points = []
            for i, j in combinations(manifest, 2):
                ki, kj = (Decimal(manifest[x]["kd_uM"]) for x in (i, j))
                if ki == kj:
                    continue
                strong, weak = (i, j) if ki < kj else (j, i)
                diff = score[arm, receptor, weak] - score[arm, receptor, strong]
                points.append(
                    dict(
                        arm=arm,
                        receptor=receptor,
                        strong=strong,
                        weak=weak,
                        measured_ddG_over_RT=abs(float((ki / kj).ln())),
                        gap=float(diff),
                        credit=float(diff > 0) + 0.5 * (diff == 0),
                        reliable=max(ki, kj) >= 3 * min(ki, kj),
                        matched=all(
                            manifest[i][k] == manifest[j][k]
                            for k in ("heavy", "charge")
                        ),
                    )
                )
            pairs.extend(points)
            kd = [float(manifest[i]["kd_uM"]) for i in manifest]
            vals = [float(score[arm, receptor, i]) for i in manifest]
            for subset, group in (
                ("all", points),
                ("3fold", [p for p in points if p["reliable"]]),
                (
                    "3fold_matched",
                    [p for p in points if p["reliable"] and p["matched"]],
                ),
            ):
                summary.append(
                    dict(
                        arm=arm,
                        receptor=receptor,
                        subset=subset,
                        pairs=len(group),
                        credit=sum(p["credit"] for p in group),
                        spearman=spearmanr(kd, vals).statistic
                        if subset == "all"
                        else "",
                        kendall_b=kendalltau(kd, vals).statistic
                        if subset == "all"
                        else "",
                    )
                )
    paired.write(out / "summary.csv", summary)
    paired.write(out / "pairs.csv", pairs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["run", "analyze"])
    parser.add_argument("out", type=Path)
    parser.add_argument("--inputs", type=Path)
    parser.add_argument(
        "--confirm",
        action="store_true",
        help="Three scientific arms, without upstream or reduced arms",
    )
    parser.add_argument("--seed", type=int, default=214039)
    parser.add_argument("--cap", type=int, default=3000)
    parser.add_argument("--cpus", nargs="+", type=int)
    for arm in ("parent", "always", "adaptive", "upstream", "counts"):
        parser.add_argument("--" + arm, type=Path)
    args = parser.parse_args()
    if args.action == "run":
        (args.out / "runs").mkdir(parents=True, exist_ok=True)
        entries = read(BENCH / "pair_gradient/galectin_external/manifest.csv")
        manifest = args.out / "settings.json"
        current = settings(args, entries)
        if manifest.is_file():
            if json.loads(manifest.read_text()) != current:
                raise ValueError("This directory contains a different experiment")
        elif any((args.out / "runs").iterdir()):
            raise ValueError("Use a fresh output directory for a new experiment")
        else:
            manifest.write_text(json.dumps(current, indent=2) + "\n")
        jobs = [(r, e["compound"]) for r in RECEPTORS for e in entries]
        random.Random(317447).shuffle(jobs)
        with ProcessPoolExecutor(len(args.cpus)) as pool:
            futures = [
                pool.submit(lane, args, jobs[i :: len(args.cpus)], cpu)
                for i, cpu in enumerate(args.cpus)
            ]
            paired.write(
                args.out / "cells.csv", [r for f in futures for r in f.result()]
            )
    summarize(args.out, args.confirm)
