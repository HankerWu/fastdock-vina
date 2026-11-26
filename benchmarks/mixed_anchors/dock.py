"""Run predeclared half-anchored trajectories and verify against pure trajectory logs."""

import argparse
import csv
import itertools
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SEEDS = (271828, 314159, 161803, 141421)


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def cases(base):
    work = base / "fdv-ionpair"
    ligands = [
        r["ligand"]
        for r in read(ROOT / "benchmarks/affinity_partition/fxa_series.csv")
        if r["included"] == "1"
    ]
    for panel, rs, ls, seeds, anchor in (
        ("cross", ("1FJS", "2P16"), ("DX9",), SEEDS, "shell"),
        (
            "trypsin",
            ("3PTB", "1BTY"),
            [r["ligand"] for r in read(work / "trypsin/series.csv")],
            SEEDS,
            "shell",
        ),
        ("inverted", ("1FJS", "2P16"), ligands, SEEDS[:2], "inverted"),
        ("fxa", ("1FJS", "2P16"), ligands, SEEDS, "shell"),
    ):
        for r, l, s in itertools.product(rs, ls, seeds):
            folder = (
                (base / ("fdv-fxa1" if r == "1FJS" else "fdv-fxa"))
                if panel != "trypsin"
                else work / "trypsin"
            )
            ligand = (work / "cross" if panel == "cross" else folder) / (l + ".pdbqt")
            config = folder / (r + ".box" if panel == "trypsin" else "box.txt")
            site = work / panel / f"{r}_{l}.{anchor}"
            prior = work / panel / f"{r}_{l}_{s}_{anchor}.log"
            native = (
                work / ("fxa" if panel == "inverted" else panel)
            ) / f"{r}_{l}_{s}_native.log"
            if panel in ("fxa", "inverted") and s in SEEDS[:2]:
                native = folder / f"{l}_{s}.log"
            yield (
                panel,
                r,
                l,
                s,
                ligand,
                folder / (r + ".pdbqt"),
                config,
                site,
                native,
                prior,
            )


def trajectories(log):
    return re.findall(r"^Search trajectory: .*", log.read_text(), re.M)


def run(case, binary, out, parity):
    panel, r, l, s, lig, rec, box, site, native, prior = case
    folder = out / panel
    folder.mkdir(exist_ok=True, parents=True)
    stem = f"{r}_{l}_{s}"
    cmd = [
        binary,
        "--ligand",
        lig,
        "--receptor",
        rec,
        "--config",
        box,
        "--seed",
        str(s),
        "--cpu",
        "1",
        "--exhaustiveness",
        "4",
        "--max_evals",
        "50000",
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--pose_selection",
        "late",
        "--search_counters",
    ]

    def run_dock(new, flags):
        with new.with_suffix(".log").open("w") as stream:
            subprocess.run(
                [*cmd, *flags, "--out", new],
                check=True,
                stdout=stream,
                stderr=subprocess.STDOUT,
            )

    if not trajectories(native):
        # First-generation native affinity outputs omit counters. Store a fresh replay in a shared
        # directory before mixed runs, so the false-head and real-head arms use the same baseline.
        replay = out / "native" / (stem + "_native.pdbqt")
        replay.parent.mkdir(exist_ok=True)
        if (
            not replay.is_file()
            or not replay.with_suffix(".log").is_file()
            or len(trajectories(replay.with_suffix(".log"))) != 4
        ):
            run_dock(replay, [])
        native = replay.with_suffix(".log")
    pool = folder / (stem + "_mixed.pdbqt")
    run_dock(pool, ["--search_anchor", site, "--search_anchor_trajectories", parity])
    n, g, m = map(trajectories, (native, prior, pool.with_suffix(".log")))
    if (
        len(n) != 4
        or len(g) != 4
        or m
        != ([n[0], g[1], n[2], g[3]] if parity == "even" else [g[0], n[1], g[2], n[3]])
    ):
        raise ValueError(f"Trajectory routing mismatch: {pool}")
    evals = sum(
        sum(map(int, re.findall(r"(?:hunt|second)_evals=(\d+)", row))) for row in m
    )
    return dict(
        panel=panel,
        receptor=r,
        ligand=l,
        seed=s,
        trajectories=len(m),
        search_evaluations=evals,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for n in ("base", "binary", "out"):
        p.add_argument(n, type=Path)
    p.add_argument("--parity", choices=("odd", "even"), default="even")
    p.add_argument("--jobs", type=int, default=4)
    a = p.parse_args()
    results = []
    # Named phases mean that native replays are complete before another prior shares one.
    for panel, group in itertools.groupby(cases(a.base), key=lambda x: x[0]):
        with ThreadPoolExecutor(a.jobs) as pool:
            results.extend(
                pool.map(lambda c: run(c, a.binary.resolve(), a.out, a.parity), group)
            )
        print(panel, sum(r["panel"] == panel for r in results), flush=True)
    write(a.out / "routing.csv", results)
