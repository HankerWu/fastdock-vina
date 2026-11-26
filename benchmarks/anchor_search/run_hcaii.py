"""Paired unrestricted, spherical-shell and directional searches on a fixed hCA II series."""

import argparse
import itertools
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "cutoff_switch" / "hcaii"))
from coordination_probe import (
    candidate,
    error_components,
    load_anchor,
    read_models,
    read_score_header,
    receptor_site,
    write_csv,
)

RECEPTORS = ("6WQ4", "6WQ8")
LIGANDS = ("6WQ4", "6WQ5", "6WQ7", "6WQ8", "6WQ9")
SEEDS = (314159, 271828)
ARMS = ("native", "shell", "point")


def dock(binary, prepared, output, cell, scaffold=None):
    rec, lig, seed, arm = cell
    stem = output / f"{rec}_{lig}_{seed}_{arm}"
    cmd = [
        str(binary),
        "--receptor",
        str(prepared / f"{rec}_receptor.pdbqt"),
        "--ligand",
        str(prepared / f"{lig}_ligand.pdbqt"),
        "--config",
        str(prepared / f"{rec}_box.txt"),
        "--cpu",
        "1",
        "--exhaustiveness",
        "4",
        "--max_evals",
        "50000",
        "--seed",
        str(seed),
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--pose_selection",
        "late",
        "--search_counters",
        "--out",
        str(stem.with_suffix(".pdbqt")),
    ]
    if arm == "scaffold":
        anchor = stem.with_suffix(".anchor")
        anchor.write_text(scaffold)
        cmd += ["--search_anchor", str(anchor)]
    elif arm != "native":
        serial, point = load_anchor(prepared, rec, lig)
        center = receptor_site(prepared, rec)["zinc_xyz"] if arm == "shell" else point
        suffix = " shell 2" if arm == "shell" else ""
        anchor = stem.with_suffix(".anchor")
        anchor.write_text(
            f"{serial} " + " ".join(f"{v:.8f}" for v in center) + suffix + "\n"
        )
        cmd += ["--search_anchor", str(anchor)]
    with stem.with_suffix(".log").open("w") as handle:
        subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=True)


def evaluate(prepared, output, cells):
    metadata = json.loads((prepared / "metadata.json").read_text())
    modes, selected, pairs = [], [], []
    for rec, lig, seed, arm in cells:
        prefix = {"rec": rec, "lig": lig, "seed": seed, "arm": arm}
        stem = output / f"{rec}_{lig}_{seed}_{arm}"
        log = stem.with_suffix(".log").read_text()
        grids = sum(
            int(x) for x in re.findall(r"(?:hunt_evals|second_evals)=(\d+)", log)
        )
        group = []
        for model, lines in read_models(
            stem.with_suffix(".pdbqt").read_text().splitlines()
        ).items():
            score = read_score_header(model, lines)["native_score"]
            atoms, info = candidate(prepared, rec, lig, lines, score, model)
            rms = error_components(
                prepared, lig, atoms, metadata[lig]["references"][rec]["coords"]
            )
            row = dict(**prefix, **info, **rms)
            modes.append(row)
            group.append(row)
        first = group[0]
        assert first["model"] == 1 and grids > 0
        for selection in ("native", "directional"):
            eligible = [r for r in group if r[selection] is not None]
            top = (
                min(eligible, key=lambda r: (r[selection], r["model"]))
                if eligible
                else None
            )
            selected.append(
                dict(
                    **prefix,
                    selection=selection,
                    grid_evals=grids,
                    best_rmsd=min(r["rmsd"] for r in group),
                    n_modes=len(group),
                    ki_nM=metadata[lig]["Ki_nM"],
                    **{
                        k: top[k] if top else None
                        for k in (
                            "model",
                            "rmsd",
                            "head_rmsd",
                            "tail_rmsd",
                            "native",
                            "zinc_N",
                            "hard",
                            "Thr199_N_to_sulfonyl_O",
                            "directional",
                        )
                    },
                )
            )
    arms = sorted({r["arm"] for r in selected})
    for rec, seed, arm, sel in itertools.product(
        RECEPTORS, SEEDS, arms, ("native", "directional")
    ):
        block = [
            r
            for r in selected
            if (r["rec"], r["seed"], r["arm"], r["selection"]) == (rec, seed, arm, sel)
        ]
        for a, b in itertools.combinations(sorted(block, key=lambda r: r["ki_nM"]), 2):
            if b["ki_nM"] < 3 * a["ki_nM"]:
                continue
            x, y = a[sel], b[sel]
            pairs.append(
                {
                    "rec": rec,
                    "seed": seed,
                    "arm": arm,
                    "selection": sel,
                    "strong": a["lig"],
                    "weak": b["lig"],
                    "score_weak_minus_strong": y - x
                    if x is not None and y is not None
                    else None,
                    "credit": (int(y > x) + 0.5 * (y == x))
                    if x is not None and y is not None
                    else None,
                }
            )
    for name, rows in [("modes", modes), ("cells", selected), ("pairs", pairs)]:
        write_csv(output / f"{name}.csv", rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("prepared", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--evaluate-only", action="store_true")
    args = parser.parse_args()
    cells = list(itertools.product(RECEPTORS, LIGANDS, SEEDS, ARMS))
    if not args.evaluate_only:
        args.output.mkdir(parents=True, exist_ok=False)
        with ThreadPoolExecutor(args.jobs) as pool:
            list(
                pool.map(
                    lambda cell: dock(
                        args.binary.resolve(), args.prepared, args.output, cell
                    ),
                    cells,
                )
            )
    evaluate(args.prepared, args.output, cells)


if __name__ == "__main__":
    main()
