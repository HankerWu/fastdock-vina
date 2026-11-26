"""Post-score check of external water maps against the withheld ligand and its observed solvent."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from evaluate import RECS, write
from scipy.spatial.distance import cdist

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from extract_independent_receptor_evidence import (
    align_donor,
    aligned_waters,
    coordinate,
    dominant_atoms,
    structure_chain,
)
from water_occupancy import box

CCD = {
    "4OGI": "R78",
    "4OGJ": "2TA",
    "4O76": "1M3",
    "4O77": "2RE",
    "6WQ4": "U7A",
    "6WQ5": "U7G",
    "6WQ7": "U6V",
    "6WQ8": "U77",
    "6WQ9": "U7J",
}


def main(a):
    records = []
    for lig, component in CCD.items():
        target = "brd4" if lig[0] == "4" else "ca2"
        work = a.brd if target == "brd4" else a.ca
        chain = "B" if lig == "4OGJ" else "A"
        model, protein = structure_chain(work / "raw" / f"{lig}.pdb", lig, chain)
        (atoms,) = [
            [
                at
                for at in dominant_atoms(res).values()
                if at.element.name not in ("D", "H")
            ]
            for ch in model
            if ch.name == chain
            for res in ch
            if res.name == component and (target == "brd4" or res.seqid.num == 303)
        ]
        pos = np.array([coordinate(at) for at in atoms])
        for rec in RECS[target]:
            dst = "B" if rec == "4OGJ" else "A"
            _, ref = structure_chain(work / "raw" / f"{rec}.pdb", rec, dst)
            center, size = box(work / "prepared" / f"{rec}_box.txt")
            _, rot, trans = align_donor(ref, protein, center, {})
            obs = np.array(
                [
                    w["xyz"]
                    for w in aligned_waters(model, lig, center, size, rot, trans)
                    if w["chain"] == chain
                ]
            )
            sites = json.loads((a.out / f"{rec}_sites.json").read_text())["sites"]
            pts = np.array([s["xyz"] for s in sites])
            ligdist = cdist(pts, pos @ rot + trans)
            nearest = ligdist.argmin(axis=1)
            d = ligdist.min(axis=1)
            w = cdist(pts, obs).min(axis=1)
            for i, s in enumerate(sites):
                records.append(
                    {
                        "target": target,
                        "query": lig,
                        "rec": rec,
                        "site_index": i,
                        "count": s["count"],
                        "closest_element": atoms[nearest[i]].element.name,
                        "closest_atom": atoms[nearest[i]].name,
                        "min_ligand": d[i],
                        "min_query_water": w[i],
                        "observed_within_1A": int(w[i] <= 1.0),
                        "crystal_displaces": int(d[i] < 2.3),
                        "within_4A": int(d[i] <= 4.0),
                    }
                )
    write(a.out / "query_site_audit.csv", records)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("brd", "ca", "out"):
        p.add_argument("--" + k, required=True, type=Path)
    main(p.parse_args())
