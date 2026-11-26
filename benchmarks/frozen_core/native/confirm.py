"""Two fresh seeds, ordinary native search/refinement, paired annotation off/on."""

import argparse
import hashlib
import itertools
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from statistics import median

from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "affinity"))
from test_entropy import atom, atom_lines, heavy, ligand_graph, pose, read_csv, write

sys.path.insert(0, str(HERE.parents[1] / "geometry_comparison" / "affinity"))
from run import STATES

SEEDS = (96431, 96457)
PREFIX = "REMARK VINA CONTACT ROTOR SCORE:"


def collect(args):
    panel = read_csv(args.inputs / "panel.csv")
    prepared = args.inputs / "prepared"
    data = args.output / "pools"
    data.mkdir(parents=True, exist_ok=True)
    rot = {
        r["compound"]: float(r["num_tors"])
        for r in read_csv(HERE.parent / "affinity" / "native_rotors.csv")
    }
    receptors = {
        p: cKDTree(
            [
                a[2]
                for a in map(atom, atom_lines((prepared / (p + ".pdbqt")).read_text()))
                if heavy(a)
            ]
        )
        for states in STATES.values()
        for p in states
    }
    graphs = {q: ligand_graph((prepared / (q + ".pdbqt")).read_text()) for q in rot}

    def experiment(job):
        row, r, seed = job
        q = row["compound"]
        key = f"{r}_{q}_{seed}"
        enabled = data / (key + ".pdbqt")
        baseline = data / (key + "-base.pdbqt")
        if args.run:
            cmd = [
                str(args.binary.resolve()),
                "--receptor",
                str(prepared / (r + ".pdbqt")),
                "--ligand",
                str(prepared / (q + ".pdbqt")),
                "--config",
                str(prepared / (r + ".box")),
                "--cpu",
                "1",
                "--exhaustiveness",
                "2",
                "--max_evals",
                "12000",
                "--num_modes",
                "9",
                "--energy_range",
                "100",
                "--seed",
                str(seed),
                "--search_counters",
            ]
            for path, opts in [(enabled, ["--contact_rotor_score"]), (baseline, [])]:
                result = subprocess.run(
                    cmd + opts + ["--out", str(path)],
                    capture_output=True,
                    text=True,
                    check=True,
                )
                path.with_suffix(".log").write_text(result.stdout)
        raw = enabled.read_text()
        base = "".join(
            l for l in raw.splitlines(keepends=True) if not l.startswith(PREFIX)
        )
        if args.run:
            assert base == baseline.read_text(), key

        def counters(log):
            return [
                l
                for l in log.read_text().splitlines()
                if l.startswith(
                    ("Search trajectory:", "Search counters:", "Hunt chain:")
                )
            ]

        assert counters(enabled.with_suffix(".log")) == counters(
            baseline.with_suffix(".log")
        )
        blocks = [b for b in raw.split("ENDMDL") if "MODEL" in b]
        differences = []
        for i, b in enumerate(blocks):
            output = next(l for l in b.splitlines() if l.startswith(PREFIX))[
                len(PREFIX) :
            ].split()
            value, ncon, neligible = float(output[0]), int(output[1]), int(output[2])
            independently = pose(b, receptors[r], graphs[q], rot[q])
            differences.append(
                (
                    ncon - independently["constrained"],
                    neligible - independently["eligible_rotors"],
                    abs(value - independently["contact"]),
                )
            )
            if not i:
                first = {
                    "native": independently["native"],
                    "raw": independently["raw"],
                    "uniform": independently["uniform"],
                    "contact": value,
                }
        return dict(
            target=row["target"],
            assay=row["assay"],
            receptor=r,
            compound=q,
            seed=seed,
            value_nM=float(row["value_nM"]),
            modes=len(blocks),
            sha256=hashlib.sha256(base.encode()).hexdigest(),
            constrained_mismatches=sum(d[0] != 0 for d in differences),
            eligible_mismatches=sum(d[1] != 0 for d in differences),
            max_error=max(d[2] for d in differences),
            **first,
        )

    tasks = [(rw, r, k) for rw in panel for r in STATES[rw["target"]] for k in SEEDS]
    with ThreadPoolExecutor(args.jobs) as executor:
        cells = list(executor.map(experiment, tasks))
    write(args.output / "cells.csv", cells)
    by = defaultdict(list)
    for c in cells:
        by[c["target"], c["assay"], c["receptor"], c["compound"]].append(c)
    groups = defaultdict(list)
    scores = ("native", "raw", "uniform", "contact")
    for (target, assay, r, q), cc in by.items():
        groups[target, assay, r].append(
            dict(
                q=q,
                v=cc[0]["value_nM"],
                **{s: median(c[s] for c in cc) for s in scores},
            )
        )
    pairs = []
    for group, compounds in groups.items():
        for left, right in itertools.combinations(compounds, 2):
            strong, weak = sorted([left, right], key=lambda v: v["v"])
            if weak["v"] < 3 * strong["v"]:
                continue
            pairs.append(
                dict(
                    target=group[0],
                    assay=group[1],
                    receptor=group[2],
                    strong=strong["q"],
                    weak=weak["q"],
                    **{
                        s: float(weak[s] > strong[s]) + 0.5 * (weak[s] == strong[s])
                        for s in scores
                    },
                )
            )
    write(args.output / "pairs.csv", pairs)
    print(
        "pools / poses",
        len(cells),
        sum(c["modes"] for c in cells),
        "native/raw/uniform/contact of",
        len(pairs),
        [sum(p[s] for p in pairs) for s in scores],
    )
    print(
        "constrained / eligible mismatches",
        *[
            sum(c[k] for c in cells)
            for k in ("constrained_mismatches", "eligible_mismatches")
        ],
        "max score error",
        max(c["max_error"] for c in cells),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--run", action="store_true")
    p.add_argument("--binary", type=Path)
    p.add_argument("--jobs", type=int, default=3)
    collect(p.parse_args())
