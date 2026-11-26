"""Score complete search matrices and compare assignment effects on identical returned poses."""

import argparse
import csv
import itertools
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from prepare import ROOT, RECEPTORS

sys.path.insert(0, str(ROOT / "tools"))
from rank_anchor_pool import ranked_pool, read_assignments
from select_zinc_pose import parse_atoms, read_models
from coordination_probe import error_components
from check_ammonium_states import ccd_vectors, handed


def save(path, rows):
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def evaluate(panel, prepared, folder, seeds, arms):
    ca = panel == "hcaii"
    metadata = json.loads(
        (prepared / ("metadata.json" if ca else "manifest.json")).read_text()
    )
    modes, cells, pairs = [], [], []
    for rec, lig, seed, arm in itertools.product(
        RECEPTORS[panel], sorted(metadata), seeds, arms
    ):
        base = folder / f"{rec}_{lig}_{seed}_{arm}"
        text = base.with_suffix(".pdbqt").read_text()
        logs = base.with_suffix(".log").read_text()
        ev = sum(int(x) for x in re.findall(r"(?:hunt_evals|second_evals)=(\d+)", logs))
        if ev < 200000:
            raise ValueError(f"Incomplete search {base}")
        scores = {}
        for prior in ("single", "multi"):
            assignments = read_assignments(
                (folder / f"{rec}_{lig}.{prior}").read_text()
            )
            scores[prior] = {
                r["model"]: r for r in ranked_pool(text, assignments)["poses"]
            }
        meta = metadata[lig]
        if ca:
            indices, reference_sign = ccd_vectors(prepared.parent, lig, meta)
        batch = []
        for num, lines in read_models(text.splitlines()).items():
            atoms = parse_atoms(lines)
            if ca:
                rms = error_components(
                    prepared, lig, atoms, meta["references"][rec]["coords"]
                )["rmsd"]
                sign = (
                    handed([np.array(atoms[s].xyz) for s in indices]) == reference_sign
                )
            else:
                xyz = np.array([atoms[s].xyz for s in meta["heavy_serials"]])
                ref = np.asarray(meta["alignments"][rec]["reference_by_match"])
                rms = float(
                    np.sqrt(np.square(ref - xyz).sum(axis=2).mean(axis=1)).min()
                )
                sign = True
            one, many = scores["single"][num], scores["multi"][num]
            row = dict(
                receptor=rec,
                ligand=lig,
                seed=seed,
                arm=arm,
                model=num,
                rmsd=rms,
                sign_ok=int(sign),
                native_affinity=one["native_affinity"],
                raw=one["native_raw"],
                single=one["landmark_penalty"],
                multi=many["landmark_penalty"],
                conditioned_single=one["guided_pose_score"],
                conditioned_multi=many["guided_pose_score"],
            )
            modes.append(row)
            batch.append(row)
        for method in (
            "released",
            "conditioned_single",
            "conditioned_multi",
            "single",
            "multi",
        ):
            top = (
                batch[0]
                if method == "released"
                else min(batch, key=lambda r: (r[method], r["model"]))
            )
            cells.append(
                dict(
                    receptor=rec,
                    ligand=lig,
                    seed=seed,
                    arm=arm,
                    method=method,
                    model=top["model"],
                    rmsd=top["rmsd"],
                    hit=int(top["rmsd"] <= 2),
                    strict=int(top["rmsd"] <= 2 and top["sign_ok"]),
                    best_rmsd=min(r["rmsd"] for r in batch),
                    strict_oracle=int(
                        any(r["rmsd"] <= 2 and r["sign_ok"] for r in batch)
                    ),
                    native_affinity=top["native_affinity"],
                    grid_evals=ev,
                )
            )
    blocks = defaultdict(list)
    for r in cells:
        if r["method"] == "released":
            blocks[(r["receptor"], r["seed"], r["arm"])].append(r)
    for (rec, seed, arm), group in blocks.items():
        get_label = lambda r: metadata[r["ligand"]]["Ki_nM" if ca else "kd_uM"]
        for a, b in itertools.combinations(sorted(group, key=get_label), 2):
            if get_label(b) < 3 * get_label(a):
                continue
            gap = b["native_affinity"] - a["native_affinity"]
            pairs.append(
                dict(
                    receptor=rec,
                    seed=seed,
                    arm=arm,
                    strong=a["ligand"],
                    weak=b["ligand"],
                    gap=gap,
                    credit=int(gap > 0) + 0.5 * int(gap == 0),
                )
            )
    for n, rs in (("modes", modes), ("cells", cells), ("pairs", pairs)):
        save(folder / f"{n}.csv", rs)
    for arm in arms:
        for method in (
            "released",
            "conditioned_single",
            "conditioned_multi",
            "single",
            "multi",
        ):
            selected = [r for r in cells if r["arm"] == arm and r["method"] == method]
            print(
                panel,
                arm,
                method,
                sum(r["hit"] for r in selected),
                sum(r["strict"] for r in selected),
                sum(r["best_rmsd"] <= 2 for r in selected),
                sum(r["strict_oracle"] for r in selected),
                round(np.median([r["rmsd"] for r in selected]), 3),
                sum(r["credit"] for r in pairs if r["arm"] == arm),
            )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("panel", choices=RECEPTORS)
    p.add_argument("prepared", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--seeds", type=int, nargs="+", default=[314159, 271828])
    p.add_argument(
        "--arms",
        nargs="+",
        choices=("native", "single", "multi"),
        default=["single", "multi"],
    )
    a = p.parse_args()
    evaluate(a.panel, a.prepared, a.out, a.seeds, a.arms)
