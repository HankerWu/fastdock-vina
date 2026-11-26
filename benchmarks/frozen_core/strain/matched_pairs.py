"""Exact equal-heavy-atom and formal-charge control, keeping receptor/assay strata."""

import argparse
from collections import defaultdict
from pathlib import Path

from pose import GEOM, read, write
from rdkit import Chem


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("affinity_work", type=Path)
    p.add_argument("native_work", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    meta = {}
    for r in read(GEOM / "affinity/panel.csv"):
        mol = Chem.MolFromSmiles(r["smiles"])
        meta[r["assay"], r["compound"]] = (
            mol.GetNumHeavyAtoms(),
            Chem.GetFormalCharge(mol),
        )
    keep = []
    for root, file, default in [
        (a.affinity_work, "control_pairs.csv", "all"),
        (a.native_work, "mmff_pairs.csv", ""),
    ]:
        for r in read(root / file):
            x = meta[r["assay"], r["stronger"]]
            if x != meta[r["assay"], r["weaker"]]:
                continue
            keep.append(
                {
                    **r,
                    "cap": r.get("cap", "fresh"),
                    "policy": r.get("policy", default),
                    "heavy": x[0],
                    "charge": x[1],
                }
            )
    # Keep a single column order for the two panels.
    first = keep[0]
    write(a.out / "matched_pair_details.csv", [{k: r[k] for k in first} for r in keep])
    group = defaultdict(list)
    for r in keep:
        for key in [(r["target"], r["assay"], r["receptor"]), ("ALL", "ALL", "ALL")]:
            group[r["cap"], *key, r["policy"], r["arm"]].append(float(r["credit"]))
    write(
        a.out / "matched_pair_summary.csv",
        [
            {
                "cap": k[0],
                "target": k[1],
                "assay": k[2],
                "receptor": k[3],
                "policy": k[4],
                "arm": k[5],
                "pairs": len(v),
                "credit": sum(v),
            }
            for k, v in sorted(group.items())
        ],
    )


if __name__ == "__main__":
    main()
