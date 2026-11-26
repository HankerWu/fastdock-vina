"""Replay the earlier two-state hCA II experiment on its exact archived coordinates."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path
import itertools
import numpy as np
from evaluate import write
from amide_rotamers import amides, dihedral_degrees
from rank_amide_pool import rank
from select_zinc_pose import read_models, parse_atoms


def compare(base, out):
    groups = defaultdict(list)
    runs = base / "fdv-map-states-declared"
    with (base / "fdv-window-fine/hcaii.csv").open() as source:
        for r in csv.DictReader(source):
            groups[r["receptor"], r["ligand"], r["seed"], r["slot"]].append(r)
    cells = defaultdict(list)
    for (rec, lig, seed, slot), g in groups.items():
        source = (
            base / "fdv-second/prepared" / f"{lig}_ligand.pdbqt"
            if slot == "common"
            else runs / "inputs" / f"{lig}_reflected.pdbqt"
        ).read_text()
        (site,) = amides(source)
        raw = (runs / f"{rec}_{lig}_{seed}_{slot}.pdbqt").read_text()
        # The one-file CLI uses printed scores; the union below recomputes its common reference.
        checked = {r["model"]: r for r in rank(source, raw, 15)}
        models = read_models(raw.splitlines())
        for r in g:
            k = int(r["model"])
            angle = dihedral_degrees(
                site, {i: a.xyz for i, a in parse_atoms(models[k]).items()}
            )
            twist = np.sin(np.radians(angle)) ** 2
            if abs(15 * twist - checked[k]["amide_penalty"]) > 1e-10:
                raise ValueError("CLI mismatch")
            cells[rec, lig, seed].append(
                dict(
                    slot=slot,
                    model=k,
                    twist=twist,
                    angle=angle,
                    inter=float(r["inter"]),
                    intra=float(r["intra"]),
                    slope=float(r["slope"]),
                    label=float(r["measurement"]),
                    rms=float(r["rmsd"]),
                    hand=int(r["sign_ok"]),
                )
            )
    rows = []
    for (rec, lig, seed), g in cells.items():
        u = min(g, key=lambda r: r["inter"] + r["intra"])["intra"]
        for k in (0, 5, 15):

            def score(r):
                return (r["inter"] + r["intra"] - u) * r["slope"] + k * r["twist"]

            t = min(g, key=score)
            rows.append(
                dict(
                    receptor=rec,
                    ligand=lig,
                    seed=seed,
                    barrier=k,
                    model=t["model"],
                    slot=t["slot"],
                    phi=t["angle"],
                    score=score(t),
                    measurement=t["label"],
                    rmsd=t["rms"],
                    strict=int(t["rms"] <= 2 and t["hand"]),
                    oracle=int(any(r["rms"] <= 2 and r["hand"] for r in g)),
                )
            )
    write(out / "historical_cells.csv", rows)
    summary = []
    for k in (0, 5, 15):
        for scope in ("all", "cross"):
            gs = [
                r
                for r in rows
                if r["barrier"] == k
                and (scope == "all" or r["receptor"] != r["ligand"])
            ]
            prs = defaultdict(list)
            for block in {(r["receptor"], r["seed"]) for r in gs}:
                batch = sorted(
                    (r for r in gs if (r["receptor"], r["seed"]) == block),
                    key=lambda r: r["measurement"],
                )
                for s, w in itertools.combinations(batch, 2):
                    if w["measurement"] >= 3 * s["measurement"]:
                        prs[s["ligand"], w["ligand"]].append(w["score"] - s["score"])
            summary.append(
                dict(
                    scope=scope,
                    barrier=k,
                    cells=len(gs),
                    top=sum(r["strict"] for r in gs),
                    oracle=sum(r["oracle"] for r in gs),
                    unique=sum((np.sign(np.median(v)) + 1) / 2 for v in prs.values()),
                    pairs=len(prs),
                )
            )
    write(out / "historical_summary.csv", summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("base", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    compare(a.base, a.out)
