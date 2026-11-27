"""Run all Gal-1/3 experimental Kd comparisons on independently generated poses."""

import argparse
import csv
import hashlib
import itertools
import random
import re
import subprocess
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from pathlib import Path
from statistics import mean

from scipy.stats import kendalltau, spearmanr

SOURCE = Path(__file__).resolve().parents[1] / "galectin_affinity"
RECEIVERS = {"Gal1": ("1GZW", "4Q26"), "Gal3": ("1KJL", "3T1L")}
SEEDS = (75467, 98251)


def save(path, rows):
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(args):
    work = args.out.resolve()
    work.mkdir(parents=True, exist_ok=True)
    with tarfile.open(SOURCE / "evidence.tar.xz") as archive:
        wanted = [
            m
            for m in archive
            if m.name == "panel.csv"
            or (m.name.startswith("prepared/") and m.name.endswith((".pdbqt", ".box")))
        ]
        archive.extractall(work, members=wanted, filter="data")
    with (work / "panel.csv").open() as file:
        panel = {int(r["compound"]): r for r in csv.DictReader(file)}
    compounds = list(range(5, 18))
    assert sorted(panel) == compounds
    binaries = {
        k: getattr(args, k).resolve() for k in ("upstream", "parent", "linear", "cubic")
    }
    if args.shift:
        binaries["shift"] = args.shift.resolve()
    save(
        work / "binaries.csv",
        [
            {"arm": k, "sha256": hashlib.sha256(v.read_bytes()).hexdigest()}
            for k, v in binaries.items()
        ],
    )
    cells = list(
        itertools.product(sum(RECEIVERS.values(), ()), compounds, SEEDS, binaries)
    )
    random.Random(872831).shuffle(cells)
    raw = work / "raw"
    raw.mkdir(exist_ok=True)
    fixtures = work / "prepared"

    def dock(cell):
        receptor, compound, seed, arm = cell
        stem = f"{receptor}-{compound}-{seed}-{arm}"
        out = raw / (stem + ".pdbqt")
        cmd = [
            binaries[arm],
            "--receptor",
            fixtures / (receptor + ".pdbqt"),
            "--ligand",
            fixtures / f"{compound}.pdbqt",
            "--config",
            fixtures / (receptor + ".box"),
            "--out",
            out,
            "--seed",
            str(seed),
            "--cpu",
            "1",
            "--max_evals",
            "12000",
            "--exhaustiveness",
            "1",
            "--num_modes",
            "9",
            "--energy_range",
            "100",
            "--verbosity",
            "0",
        ]
        if arm != "upstream":
            cmd += ["--search_counters"]
        start = time.perf_counter()
        log = subprocess.run(cmd, capture_output=True, check=True, text=True).stdout
        seconds = time.perf_counter() - start
        (raw / (stem + ".log")).write_text(log)
        counters = {
            k: int(v)
            for k, v in re.findall(
                r"\b(hunt_evals|second_evals|explicit_evals|hunt_calls|second_calls|explicit_calls)=(\d+)",
                log,
            )
        }
        if len(counters) != 6:
            raise ValueError(
                "Run the counts-instrumented upstream build for this panel: " + stem
            )
        energies = re.findall(
            r"^REMARK VINA RESULT:\s+([-\d.]+)", out.read_text(), re.MULTILINE
        )
        return dict(
            receptor=receptor,
            compound=compound,
            seed=seed,
            arm=arm,
            modes=len(energies),
            top_score=float(energies[0]),
            seconds=round(seconds, 6),
            **counters,
        )

    if args.extend_shift:
        if not args.shift:
            raise ValueError("--extend-shift requires --shift")
        with (work / "cells.csv").open() as file:
            rows = list(csv.DictReader(file))
        expected = {c for c in cells if c[3] != "shift"}
        obtained = {
            (r["receptor"], int(r["compound"]), int(r["seed"]), r["arm"]) for r in rows
        }
        if obtained != expected or len(rows) != len(expected):
            raise ValueError("Incomplete or previously extended four-arm panel")
        for r in rows:
            r.update(
                compound=int(r["compound"]),
                seed=int(r["seed"]),
                top_score=float(r["top_score"]),
            )
        cells = [c for c in cells if c[3] == "shift"]
    else:
        rows = []
    with ThreadPoolExecutor(args.jobs) as pool:
        rows.extend(pool.map(dock, cells))
    save(work / "cells.csv", rows)
    indexed = {
        (r["receptor"], r["compound"], r["seed"], r["arm"]): r["top_score"]
        for r in rows
    }
    pairs, summary = [], []
    for target, states in RECEIVERS.items():
        kd = {i: Decimal(panel[i][target + "_Kd_uM"]) for i in compounds}
        for arm in binaries:
            for label in (*states, "consensus"):
                used = states if label == "consensus" else (label,)
                scores = {
                    i: mean(
                        indexed[state, i, seed, arm] for state in used for seed in SEEDS
                    )
                    for i in compounds
                }
                for subset in ("all", "triazoles", "equal_heavy"):
                    use = compounds[:-1] if subset == "triazoles" else compounds
                    correct = []
                    for a, b in itertools.combinations(use, 2):
                        a, b = sorted((a, b), key=kd.get)
                        if kd[b] < 3 * kd[a]:
                            continue
                        if (
                            subset == "equal_heavy"
                            and panel[a]["heavy_atoms"] != panel[b]["heavy_atoms"]
                        ):
                            continue
                        credit = (
                            1
                            if scores[a] < scores[b]
                            else (0.5 if scores[a] == scores[b] else 0)
                        )
                        correct.append(credit)
                        pairs.append(
                            {
                                "target": target,
                                "state": label,
                                "arm": arm,
                                "subset": subset,
                                "strong": a,
                                "weak": b,
                                "strong_score": scores[a],
                                "weak_score": scores[b],
                                "credit": credit,
                            }
                        )
                    rho = tau = ""
                    if subset != "equal_heavy":
                        x, y = [scores[i] for i in use], [float(kd[i]) for i in use]
                        rho, tau = spearmanr(x, y).statistic, kendalltau(x, y).statistic
                    summary.append(
                        {
                            "target": target,
                            "state": label,
                            "arm": arm,
                            "subset": subset,
                            "correct": sum(correct),
                            "pairs": len(correct),
                            "spearman": rho,
                            "kendall_tau_b": tau,
                        }
                    )
    save(work / "pairs.csv", pairs)
    save(work / "summary.csv", summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("upstream", "parent", "linear", "cubic", "out"):
        parser.add_argument(name, type=Path)
    parser.add_argument(
        "--shift", type=Path, help="Follow-up explicit shifted-cutoff build"
    )
    parser.add_argument(
        "--extend-shift",
        action="store_true",
        help="Reuse complete four-arm output in out",
    )
    parser.add_argument("--jobs", type=int, default=1)
    run(parser.parse_args())
