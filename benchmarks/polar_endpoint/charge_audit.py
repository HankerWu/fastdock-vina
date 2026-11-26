"""Report within-charge and charge-changing pair credit with the prepared charge state."""

import argparse
import statistics
from collections import defaultdict
from pathlib import Path

from polar_endpoint import atoms
from run import read, write


def audit(results, locations):
    ligands = []
    charges = {}
    for target, folder in locations.items():
        scores = read(results / (target + ".csv"))
        for lig in sorted({r["ligand"] for r in scores}):
            group = [r for r in scores if r["ligand"] == lig]
            file = folder / (
                f"{lig}.pdbqt"
                if target.startswith("ca_holdout")
                else f"{lig}_ligand.pdbqt"
            )
            q = float(sum(atoms(file.read_text().splitlines()).charges))
            if abs(q - round(q)) > 0.05:
                raise ValueError(f"Ambiguous integer charge {file}: {q}")
            charges[target, lig] = round(q)
            ligands.append(
                {
                    "target": target,
                    "ligand": lig,
                    "prepared_charge": q,
                    "poses": len(group),
                    "median_coulomb": statistics.median(
                        float(r["coulomb"]) for r in group
                    ),
                    "median_polarization": statistics.median(
                        float(r["delta_gb"]) for r in group
                    ),
                    "median_transfer": statistics.median(
                        float(r["polar_transfer"]) for r in group
                    ),
                }
            )
    groups = defaultdict(list)
    for r in read(results / "pairs.csv"):
        target = r["target"]
        same = charges[target, r["strong"]] == charges[target, r["weak"]]
        groups[target, r["pool"], r["method"], same].append(r)
    rows = []
    for (target, pool, method, same), pairs in sorted(groups.items()):
        distinct = defaultdict(list)
        for p in pairs:
            distinct[p["strong"], p["weak"]].append(float(p["diff"]))
        med = [statistics.median(v) for v in distinct.values()]
        rows.append(
            {
                "target": target,
                "pool": pool,
                "method": method,
                "same_charge": int(same),
                "pairs": len(pairs),
                "correct": sum(float(r["correct"]) for r in pairs),
                "unique": len(med),
                "unique_correct": sum((x > 0) + 0.5 * (x == 0) for x in med),
            }
        )
    return ligands, rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results", type=Path)
    for field in ("hcaii", "brd4", "t4l", "ca_holdout36", "ca_holdout20"):
        p.add_argument("--" + field, required=True, type=Path)
    a = p.parse_args()
    locations = vars(a).copy()
    locations.pop("results")
    ligands, rows = audit(a.results, locations)
    write(a.results / "ligands.csv", ligands)
    write(a.results / "charge_pairs.csv", rows)
