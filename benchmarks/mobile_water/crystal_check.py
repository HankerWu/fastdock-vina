"""Post-prediction check of the five thrombin co-crystals, in their own frames."""

import argparse
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from mobile_water import PRESENT, RT, Site
from replay import pdb_lines, records, save, thrombin_site, xyz


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--throm", type=Path, required=True)
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    assay = records(args.repo / "benchmarks/thrombin_water/expected/series.csv")
    raw = args.throm / "raw"
    external = pdb_lines(raw / "3D49.pdb")
    evidence = []
    for entry in assay:
        pdb = entry["PDB"]
        if not pdb:
            continue
        r = pdb_lines(raw / f"{pdb}.pdb")
        _, center, b = thrombin_site(external, r)[0]
        protein = cKDTree([xyz(l) for l in r if l.startswith("ATOM") and l[21] in "HL"])
        ligand = cKDTree([xyz(l) for l in r if l[17:20] == entry["CCD"]])
        waters = cKDTree([xyz(l) for l in r if l[17:20] == "HOH"])
        rigid, mobile = Site.prepare(center, b, protein).evaluate(ligand)
        # If F=-RT log(Z), posterior P(occupied)=1-(1-p)/Z.
        evidence.append(
            {
                "ligand": entry["id"],
                "pdb": pdb,
                "ligand_min": ligand.query(center)[0],
                "observed_water_min": waters.query(center)[0],
                "rigid_penalty": rigid,
                "mobile_penalty": mobile,
                "rigid_posterior": 1 - (1 - PRESENT) * np.exp(rigid / RT),
                "mobile_posterior": 1 - (1 - PRESENT) * np.exp(mobile / RT),
            }
        )
    save(args.out, evidence)


if __name__ == "__main__":
    main()
