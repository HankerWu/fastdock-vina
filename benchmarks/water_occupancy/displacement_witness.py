"""Test independent published BRD4 ligands as witnesses that specific waters are displaceable."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from evaluate import RECS, brd_pools, parse_pool, read, score_components, write
from water_occupancy import box, displacement_witnesses

WITNESSES = {"6MH1": "JQP", "5I88": "69G", "4GPJ": "0Q1", "5DLZ": "5D1"}


def run(a):
    work = a.brd
    out = a.out
    negative = {}
    manifest = {}
    mask_rows = []
    for rec in RECS["brd4"]:
        sites = json.loads((out / f"{rec}_sites.json").read_text())["sites"]
        chain = "B" if rec == "4OGJ" else "A"
        center, _ = box(work / "prepared" / (rec + "_box.txt"))
        sources = [
            (name, "A", component, a.donors / (name + ".pdb"))
            for name, component in WITNESSES.items()
        ]
        for name, _, _, path in sources:
            manifest[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        supports = displacement_witnesses(
            work / "raw" / f"{rec}.pdb", rec, chain, center, sites, sources
        )
        for i, s in enumerate(sites):
            mask_rows.append(
                {
                    "rec": rec,
                    "index": i,
                    "positive_structures": s["count"],
                    "displacement_witnesses": "+".join(supports[i]),
                }
            )
        negative[rec] = [s for i, s in enumerate(sites) if supports[i]]
    write(out / "witness_sites.csv", mask_rows)
    (out / "witness_donors.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    saved = {
        (r["rec"], r["lig"], str(r["seed"]), r["arm"], r["candidate"]): r
        for r in read(out / "poses.csv")
        if r["panel"] == "brd4"
    }
    cells = []
    pairs = []
    for _, rec, lig, seed, arm, label, models in brd_pools(a.brd_pools):
        rows = []
        for mid, path, index, _, rms, _ in models:
            orig = saved[rec, lig, str(seed), arm, mid]
            delta = score_components(negative[rec], *parse_pool(path)[index])["joint"]
            rows.append(
                {"candidate": mid, "score": float(orig["joint"]) - delta, "rmsd": rms}
            )
        chosen = min(rows, key=lambda r: r["score"])
        cells.append(
            {"rec": rec, "lig": lig, "seed": seed, "arm": arm, "label": label, **chosen}
        )
    write(out / "witness_cells.csv", cells)
    for arm, group in itertools.groupby(
        sorted(cells, key=lambda r: r["arm"]), lambda r: r["arm"]
    ):
        group = list(group)
        credit = 0
        for (rec, seed), block in itertools.groupby(
            sorted(group, key=lambda r: (r["rec"], r["seed"])),
            lambda r: (r["rec"], r["seed"]),
        ):
            for x, y in itertools.combinations(
                sorted(block, key=lambda r: r["label"]), 2
            ):
                if y["label"] < 3 * x["label"]:
                    continue
                gap = y["score"] - x["score"]
                correct = int(gap > 0) + 0.5 * (gap == 0)
                credit += correct
                pairs.append(
                    {
                        "rec": rec,
                        "seed": seed,
                        "arm": arm,
                        "strong": x["lig"],
                        "weak": y["lig"],
                        "gap": gap,
                        "credit": correct,
                    }
                )
        print(arm, sum(r["rmsd"] <= 2 for r in group), credit)
    write(out / "witness_pairs.csv", pairs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("donors", "brd", "brd_pools", "out"):
        p.add_argument("--" + k.replace("_", "-"), type=Path, required=True)
    run(p.parse_args())
