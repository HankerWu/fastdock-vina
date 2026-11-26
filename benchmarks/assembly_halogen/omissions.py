"""Compare deposited holo components to the protein-only experiment."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import csv
import sys
import urllib.request

import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist


def points(lines):
    return np.array([[float(l[i : i + 8]) for i in (30, 38, 46)] for l in lines])


def audit(case, root):
    folder = root / case
    pdb, cc = case.split("_")
    raw = folder / "deposition.pdb"
    if not raw.exists():
        raw.write_bytes(
            urllib.request.urlopen(
                f"https://files.rcsb.org/download/{pdb}.pdb", timeout=45
            ).read()
        )
    x = (
        Chem.SDMolSupplier(str(folder / "crystal.sdf"), removeHs=True)[0]
        .GetConformer()
        .GetPositions()
    )
    groups = {}
    for l in raw.read_text().splitlines():
        if (
            l.startswith("HETATM")
            and l[16] in " A"
            and l[76:78].strip() not in ("H", "D")
        ):
            groups.setdefault(l[17:27].strip(), []).append(l)
    nearest = {g: cdist(x, points(ls)).min() for g, ls in groups.items()}
    matched = min((g for g in groups if g.split()[0] == cc), key=nearest.get)
    close = sorted(
        (
            (float(d), g)
            for g, d in nearest.items()
            if g != matched and d < 6.0 and g.split()[0] not in ("HOH", "DOD")
        ),
        key=lambda v: v[0],
    )
    return dict(
        case=case,
        crystal=matched,
        nearest_native_distance=round(float(nearest[matched]), 3),
        nonwater_hetero="; ".join(f"{g} ({d:.2f} A)" for d, g in close),
    )


if __name__ == "__main__":
    work = Path(sys.argv[1])
    with ThreadPoolExecutor(4) as pool:
        rows = list(
            pool.map(
                lambda case: audit(case, work),
                Path(__file__).with_name("panel.txt").read_text().split(),
            )
        )
    with open(sys.argv[2], "w") as out:
        w = csv.DictWriter(out, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)
