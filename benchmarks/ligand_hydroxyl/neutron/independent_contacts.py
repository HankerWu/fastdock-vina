"""Independently enumerate nearest periodic images about the deposited ligand O-Ds."""

import argparse
import csv
import gzip
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def crystal_contacts(code):
    lines = (
        gzip.decompress((HERE / f"{code}.pdb.gz").read_bytes()).decode().splitlines()
    )
    cell = next(l for l in lines if l.startswith("CRYST1"))
    lengths = np.array([float(cell[i:j]) for i, j in [(6, 15), (15, 24), (24, 33)]])
    if [float(cell[i:j]) for i, j in [(33, 40), (40, 47), (47, 54)]] != [
        90.0,
        90.0,
        90.0,
    ]:
        raise ValueError("This audit uses the orthorhombic cells deposited here")
    operators = {}
    for line in lines:
        if line.startswith("REMARK 290   SMTRY"):
            tag, num, *values = line.split()[2:]
            operators.setdefault(int(num), {})[int(tag[-1]) - 1] = [
                float(x) for x in values
            ]
    ops = [np.array([op[i] for i in range(3)]) for _, op in sorted(operators.items())]
    if not len(ops) == 4:
        raise ValueError("Incomplete explicit space-group operators")
    atoms = [
        l
        for l in lines
        if l.startswith(("ATOM  ", "HETATM"))
        and float(l[54:60]) > 0
        and l[76:78].strip() in ("O", "N")
    ]
    xyz = np.array([[float(l[s : s + 8]) for s in (30, 38, 46)] for l in atoms])
    donor_ligands = {
        (l[17:27], l[12:16].strip()[1:])
        for l in lines
        if l.startswith("HETATM")
        and l[76:78].strip() == "D"
        and l[17:20] in ("GOL", "BGC", "GAL")
        and l[12:16].strip().startswith("DO")
    }
    rows = []
    for k, line in enumerate(atoms):
        if (line[17:27], line[12:16].strip()) not in donor_ligands:
            continue
        alternate = [
            (np.linalg.norm(xyz[k] - xyz[j]), atoms[j][12:27].strip())
            for j in range(len(atoms))
            if atoms[j][16] not in " A"
        ]
        nearest_alt = min(alternate, default=(np.inf, ""))
        nearest = None
        for opid, op in enumerate(ops):
            transformed = xyz @ op[:, :3].T + op[:, 3]
            # The nearest integer fractional displacement suffices in an orthorhombic cell.
            shifts = np.rint((xyz[k] - transformed) / lengths).astype(int)
            pts = transformed + lengths * shifts
            d = np.linalg.norm(pts - xyz[k], axis=1)
            if opid == 0:
                d[(shifts == 0).all(axis=1)] = np.inf
            j = np.argmin(d)
            cand = (
                d[j],
                atoms[j][12:27].strip(),
                opid + 1,
                ",".join(map(str, shifts[j])),
            )
            if nearest is None or cand[0] < nearest[0]:
                nearest = cand
        # Any second-nearest image of an atom has distance at least half
        # the shortest unit-cell edge, including translated copies whose
        # nearest image above was the discarded identity.
        if nearest[0] >= lengths.min() / 2:
            raise ValueError("Expand the periodic image search for this cell")
        rows.append(
            dict(
                pdb=code,
                site=line[12:27].strip(),
                nearest_mate_A=nearest[0],
                atom=nearest[1],
                operator=nearest[2],
                unit_translations=nearest[3],
                active_within_3p5=int(nearest[0] < 3.5),
                nearest_excluded_ASU_alt_A=nearest_alt[0],
                nearest_excluded_ASU_alt_atom=nearest_alt[1],
            )
        )
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    rows = crystal_contacts("6EXY") + crystal_contacts("6EYM")
    with a.output.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
