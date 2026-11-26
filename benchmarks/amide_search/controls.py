"""Compare the amide penalty with the same penalty on its adjacent single bond."""

import argparse
import csv
import itertools
from collections import defaultdict
from pathlib import Path
import numpy as np
from rdkit import Chem
from evaluate import write
from amide_rotamers import amides, dihedral_degrees
from make_scaffold_anchors import query_graph
from select_zinc_pose import read_models, parse_atoms


def read(file):
    with file.open() as f:
        return list(csv.DictReader(f))


def check(base, folder, out):
    new = [
        r for r in read(out / "modes.csv") if r["panel"] == "ca" and r["form"] == "free"
    ]
    old = read(base / "fdv-window-fine/hcaii.csv")
    pools = defaultdict(list)
    for r in new:
        stem = f"ca_{r['rec']}_{r['ligand']}_{r['state']}_{r['seed']}_{r['slot']}_free"
        pools[
            "new", r["rec"], r["ligand"], r["seed"], folder / (stem + ".pdbqt")
        ].append(r)
    for r in old:
        stem = f"{r['receptor']}_{r['ligand']}_{r['seed']}_{r['slot']}"
        pools[
            "archive",
            r["receptor"],
            r["ligand"],
            r["seed"],
            base / "fdv-map-states-declared" / (stem + ".pdbqt"),
        ].append(r)
    cells = defaultdict(list)
    for (cohort, rec, lig, seed, path), rows in pools.items():
        src = (base / "fdv-second/prepared" / f"{lig}_ligand.pdbqt").read_text()
        (amide,) = amides(src)
        mol, serials = query_graph(src)
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[O]=[CX3]-[CH2]-[#7+]"))
        links = [
            tuple(serials[i] for i in m) for m in matches if serials[m[1]] == amide[1]
        ]
        (neighbor,) = links
        models = read_models(path.read_text().splitlines())
        for r in rows:
            key = int(r.get("model", r.get("mode")))
            xyz = {i: a.xyz for i, a in parse_atoms(models[key]).items()}
            x, y = (np.radians(dihedral_degrees(s, xyz)) for s in (amide, neighbor))
            cells[cohort, rec, lig, seed].append(
                dict(
                    raw=float(r["inter"]) + float(r["intra"]),
                    intra=float(r["intra"]),
                    slope=float(r["slope"]),
                    amide=np.sin(x) ** 2,
                    adjacent=np.sin(y) ** 2,
                    perpendicular=np.cos(x) ** 2,
                    zero=0.0,
                    rms=float(r["rmsd"]),
                    hand=int(r.get("sign_ok", r.get("hand"))),
                    value=float(r["measurement"]),
                )
            )
    tops = []
    for (cohort, rec, lig, seed), g in cells.items():
        u = min(g, key=lambda r: r["raw"])["intra"]
        for variant in ("zero", "amide", "adjacent", "perpendicular"):

            def energy(r):
                return (r["raw"] - u) * r["slope"] + 15 * r[variant]

            w = min(g, key=energy)
            tops.append(
                dict(
                    cohort=cohort,
                    rec=rec,
                    ligand=lig,
                    seed=seed,
                    variant=variant,
                    score=energy(w),
                    hit=int(w["rms"] <= 2 and w["hand"]),
                    rmsd=w["rms"],
                    measurement=w["value"],
                )
            )
    write(out / "torsion_controls.csv", tops)
    summary = []
    for cohort, scope, variant in itertools.product(
        ("new", "archive"),
        ("all", "cross"),
        ("zero", "amide", "adjacent", "perpendicular"),
    ):
        g = [
            r
            for r in tops
            if r["cohort"] == cohort
            and r["variant"] == variant
            and (scope == "all" or r["ligand"] != r["rec"])
        ]
        edges = defaultdict(list)
        for block in {(r["rec"], r["seed"]) for r in g}:
            blockrows = sorted(
                (r for r in g if (r["rec"], r["seed"]) == block),
                key=lambda r: r["measurement"],
            )
            for a, b in itertools.combinations(blockrows, 2):
                if b["measurement"] >= 3 * a["measurement"]:
                    edges[a["ligand"], b["ligand"]].append(b["score"] - a["score"])
        summary.append(
            dict(
                cohort=cohort,
                scope=scope,
                variant=variant,
                hits=sum(r["hit"] for r in g),
                cells=len(g),
                median_rmsd=np.median([r["rmsd"] for r in g]),
                unique=sum((np.sign(np.median(e)) + 1) / 2 for e in edges.values()),
                pairs=len(edges),
            )
        )
    write(out / "torsion_control_summary.csv", summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("base", "folder", "out"):
        p.add_argument(name, type=Path)
    args = p.parse_args()
    check(args.base, args.folder, args.out)
