"""Prepare both deposited MerTK chains and keep every omitted heavy-atom distance auditable."""

import argparse
import json
from pathlib import Path
import sys
import time
from urllib.request import urlopen
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("work", type=Path)
a = p.parse_args()
here = Path(__file__).resolve().parent
sys.path.insert(0, str(here.parents[2] / "benchmarks/paired_crossdocking"))
import run as base

start = time.perf_counter()
raw = a.work / "raw"
d = a.work / "5K0X"
meta = json.loads((d / "preparation.json").read_text())
lines = (raw / "5K0X.pdb").read_text().splitlines()
xyz = np.array(meta["mapped"][0])
omit = [
    l
    for l in lines
    if l.startswith(("ATOM  ", "HETATM"))
    and l[17:20].strip() != "HOH"
    and l[76:78].strip() not in ("H", "D")
    and not (l.startswith("ATOM  ") and l[21] in ("A", "B"))
    and not (l[21] == "A" and l[17:20] == "K0X")
]
closest = min(
    (float(np.linalg.norm(xyz - base.coords(l), axis=1).min()), l[12:27]) for l in omit
)
if closest[0] < 8:
    raise ValueError(closest)
reports = {}
atoms = []
for chain in ("A", "B"):
    group = [l for l in lines if l.startswith("ATOM  ") and l[21] == chain]
    numbered = {
        rid: i for i, rid in enumerate(dict.fromkeys(l[22:27] for l in group), 1)
    }
    transformed = [l[:21] + "A" + f"{numbered[l[22:27]]:4d} " + l[27:] for l in group]
    sub = d / chain
    sub.mkdir(exist_ok=True)
    report = base.receptor(transformed, sub, "5K0X", xyz)
    near = report["closest_repaired_ligand_distance"]
    if near is not None and near < 8:
        raise ValueError((chain, near))
    reports[chain] = report
    for l in (sub / "receptor.pdbqt").read_text().splitlines():
        if l.startswith(("ATOM  ", "HETATM")):
            atoms.append(l[:6] + f"{len(atoms) + 1:5d}" + l[11:21] + chain + l[22:])
(d / "receptor.pdbqt").write_text("\n".join(atoms) + "\n")
meta["assembly"] = {
    "chains": reports,
    "nearest_omitted_nonwater": closest,
    "extra_preparation_s": time.perf_counter() - start,
}
(d / "assembly.json").write_text(json.dumps(meta["assembly"], indent=2) + "\n")
sdf = raw / "5K0X-bound.sdf"
if not sdf.exists():
    sdf.write_bytes(
        urlopen(
            "https://models.rcsb.org/v1/5k0x/ligand?auth_asym_id=A&auth_seq_id=904&encoding=sdf&copy_all_categories=false",
            timeout=45,
        ).read()
    )
