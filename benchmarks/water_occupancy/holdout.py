"""Two-box replication on other zinc isoforms with newly held-out water sources."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from evaluate import (
    atom_arrays,
    box,
    build_sites,
    read,
    read_models,
    score_components,
    write,
)

DONORS = {
    "CAI": ("1HCB", "1HUG", "2CAB", "2FOY"),
    "CAXII": ("1JCZ", "4WW8", "5MSA", "4QJ0"),
}


def run(args):
    out = args.out
    out.mkdir(exist_ok=True, parents=True)
    hashes = {}
    candidates = []
    selections = []
    pairs = []
    summary = []
    for width, root in [(36, args.large), (20, args.small)]:
        prep = root / "prepared"
        key = lambda r: (r["receptor"], r["ligand"], r["crystal"], r["seed"], r["arm"])
        all_sites = {}
        for rec in ("1AZM", "1BZM", "1JD0"):
            target = "CAXII" if rec == "1JD0" else "CAI"
            names = DONORS[target]
            donors = [(n, "A", args.donors / f"{n}.pdb") for n in names]
            for name, _, p in donors:
                hashes[name] = hashlib.sha256(p.read_bytes()).hexdigest()
            info = build_sites(
                root / "raw" / (rec + ".pdb"),
                rec,
                "A",
                prep / (rec + "_receptor.pdbqt"),
                *box(prep / (rec + ".box")),
                donors,
            )
            (out / f"{target}_{rec}_{width}_sites.json").write_text(
                json.dumps(info, indent=2) + "\n"
            )
            all_sites[rec] = info["sites"]
        for (rec, lig, crystal, seed, arm), group in itertools.groupby(
            sorted(read(root / "holdout_modes.csv"), key=key), key
        ):
            actual = {int(r["model"]): r for r in group}
            f = root / "runs" / f"{rec}_{crystal}_{seed}_{arm}.pdbqt"
            q = {
                "box": width,
                "target": "CAXII" if rec == "1JD0" else "CAI",
                "receptor": rec,
                "ligand": lig,
                "seed": seed,
                "arm": arm,
            }
            rows = []
            for i, lines in read_models(f.read_text().splitlines()).items():
                xyz, pol, _ = atom_arrays(lines)
                native = float(actual[i]["score"])
                terms = score_components(all_sites[rec], xyz, pol)
                row = {
                    **q,
                    "model": i,
                    "rmsd": float(actual[i]["rmsd"]),
                    "native": native,
                    **{k: native + v for k, v in terms.items()},
                }
                candidates.append(row)
                rows.append(row)
            for method in (
                "native",
                "exclusion",
                "contact",
                "joint",
                "replace_any",
                "replace_polar",
            ):
                win = min(rows, key=lambda r: r[method])
                selections.append(
                    {
                        **q,
                        "method": method,
                        "model": win["model"],
                        "rmsd": win["rmsd"],
                        "score": win[method],
                        "oracle": min(r["rmsd"] for r in rows),
                    }
                )

    def outer(r):
        return r["box"], r["target"], r["arm"], r["method"]

    for (width, target, arm, method), group in itertools.groupby(
        sorted(selections, key=outer), key=outer
    ):
        group = list(group)
        key = lambda r: (r["receptor"], r["seed"])
        credits = []
        if target == "CAI":
            for (rec, seed), block in itertools.groupby(
                sorted(group, key=key), key=key
            ):
                vals = {r["ligand"]: r["score"] for r in block}
                d = vals["AZM"] - vals["MZM"]
                credits.append(int(d > 0) + 0.5 * (d == 0))
                pairs.append(
                    {
                        "box": width,
                        "arm": arm,
                        "method": method,
                        "rec": rec,
                        "seed": seed,
                        "gap": d,
                        "credit": credits[-1],
                    }
                )
        summary.append(
            {
                "box": width,
                "target": target,
                "arm": arm,
                "method": method,
                "n": len(group),
                "top_2A": sum(r["rmsd"] <= 2 for r in group),
                "pool_2A": sum(r["oracle"] <= 2 for r in group),
                "median": np.median([r["rmsd"] for r in group]),
                "credit": sum(credits),
                "pairs": len(credits),
            }
        )
    (out / "donors.json").write_text(
        json.dumps(hashes, sort_keys=True, indent=2) + "\n"
    )
    for n, rows in (
        ("candidates", candidates),
        ("selections", selections),
        ("pairs", pairs),
        ("summary", summary),
    ):
        write(out / (n + ".csv"), rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for a in ("donors", "large", "small", "out"):
        p.add_argument("--" + a, required=True, type=Path)
    run(p.parse_args())
