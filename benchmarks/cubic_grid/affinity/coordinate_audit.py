"""Independently compare input and prepared heavy atoms and report unresolved segment bounds."""

import json
import sys
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist

root = Path(sys.argv[1])
out = {}


def xyz(l):
    return [float(l[k : k + 8]) for k in (30, 38, 46)]


for state in ("6YT6", "6ZCY"):
    lines = (root / "raw" / f"{state}.pdb").read_text().splitlines()
    lig = np.array([xyz(l) for l in lines if l.startswith("HETATM")])
    p = [
        l for l in lines if l.startswith("ATOM") and l[76:78].strip() not in ("H", "D")
    ]
    q = [
        l
        for l in (root / "prepared" / f"{state}-protein.pdb").read_text().splitlines()
        if l.startswith("ATOM") and l[76:78].strip() not in ("H", "D")
    ]
    atoms = lambda rows: {(int(l[22:26]), l[12:16].strip()): xyz(l) for l in rows}
    old, new = atoms(p), atoms(q)
    common = old.keys() & new.keys()
    max_move = max(float(np.linalg.norm(np.subtract(old[k], new[k]))) for k in common)
    added = sorted(new.keys() - old.keys())
    d = float(cdist([new[k] for k in added], lig).min()) if added else None
    cas = [k[0] for k in old if k[1] == "CA"]
    gaps = [
        {
            "after": a,
            "before": b,
            "nearest_flanking_CA_A": float(
                cdist([old[a, "CA"], old[b, "CA"]], lig).min()
            ),
        }
        for a, b in zip(cas, cas[1:])
        if b > a + 1
    ]
    severe = []
    for k in added:
        exclude = [j for j in new if abs(j[0] - k[0]) > 1]
        distances = cdist([new[k]], [new[j] for j in exclude])[0]
        j = int(distances.argmin())
        if distances[j] < 2:
            severe.append(
                {"added": k, "other": exclude[j], "distance_A": float(distances[j])}
            )
    out[state] = {
        "shared_atom_max_displacement_A": max_move,
        "added_heavy_atom_count": len(added),
        "nearest_added_atom_A": d,
        "residue_number_gaps": gaps,
        "contacts_below_2A_nonadjacent_residues": severe,
    }
(root / "coordinate_audit.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out))
