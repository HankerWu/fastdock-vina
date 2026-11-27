"""First-pose side-chain displacement and chi1 rotation in the rigid receptor frame."""

import argparse
import csv
from pathlib import Path

import numpy as np

# Symmetries of the prepared states. Paired aromatic swaps preserve connectivity.
SWAPS = {
    "GLU": (("OE1", "OE2"),),
    "ASP": (("OD1", "OD2"),),
    "ARG": (("NH1", "NH2"),),
    "TYR": (("CD1", "CD2"), ("CE1", "CE2")),
    "PHE": (("CD1", "CD2"), ("CE1", "CE2")),
}


def coordinates(line):
    return np.array([float(line[i : i + 8]) for i in (30, 38, 46)])


def read_sidechain(text):
    first = text.split("ENDMDL", 1)[0]
    if first.count("BEGIN_RES ") != 1:
        raise ValueError("The panel specifies one flexible residue")
    domain = first.split("BEGIN_RES ", 1)[1].split("END_RES ", 1)[0].splitlines()
    identifier = domain[0].split()
    xyz = {}
    for line in domain:
        if line.startswith(("ATOM", "HETATM")) and line[77:].strip() not in {"HD", "H"}:
            name = line[12:16].strip()
            if name in xyz:
                raise ValueError(f"Duplicate atom {name}")
            xyz[name] = coordinates(line)
    return identifier, xyz


def chi1(backbone_n, coords):
    fourth = next(k for k in ("OG1", "OG", "CG") if k in coords)
    a = backbone_n - coords["CA"]
    axis = coords["CB"] - coords["CA"]
    axis /= np.linalg.norm(axis)
    b = coords[fourth] - coords["CB"]
    a -= a.dot(axis) * axis
    b -= b.dot(axis) * axis
    return float(np.rad2deg(np.arctan2(np.cross(axis, a).dot(b), a.dot(b))))


def measure(prepared, mode):
    identifier, before = read_sidechain((prepared / "partition_flex.pdbqt").read_text())
    actual_id, after = read_sidechain(mode.read_text())
    if identifier != actual_id or before.keys() != after.keys():
        raise ValueError(f"Inconsistent flexible residue for {mode}")
    if any(np.linalg.norm(before[k] - after[k]) > 0.002 for k in ("CA", "CB")):
        raise ValueError(f"Fixed CA-CB axis moved in {mode}")
    mobile = sorted(before.keys() - {"CA", "CB"})
    mapping = dict(zip(mobile, mobile))

    def rmsd():
        return float(
            np.sqrt(
                np.mean([np.sum((before[k] - after[mapping[k]]) ** 2) for k in mobile])
            )
        )

    direct = rmsd()
    for a, b in SWAPS.get(identifier[0], ()):
        mapping[a], mapping[b] = b, a
    alternate = rmsd()
    n_atoms = [
        line
        for line in (prepared / "partition_rigid.pdbqt").read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
        and line[12:16].strip() == "N"
        and [line[17:20].strip(), line[21:22].strip(), line[22:26].strip()]
        == identifier
    ]
    if len(n_atoms) != 1:
        raise ValueError(f"Cannot identify backbone N for {prepared}")
    n = coordinates(n_atoms[0])
    initial, final = chi1(n, before), chi1(n, after)
    return {
        "residue": " ".join(identifier),
        "mobile_heavy_atoms": len(mobile),
        "mobile_atom_names": "/".join(mobile),
        "top1_mobile_sc_rmsd": min(direct, alternate),
        "top1_name_mapped_rmsd": direct,
        "symmetric_swap_used": alternate < direct,
        "input_chi1_deg": initial,
        "output_chi1_deg": final,
        "absolute_chi1_displacement_deg": abs((final - initial + 180) % 360 - 180),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("output", type=Path)
    ap.add_argument(
        "--panel",
        nargs=2,
        type=Path,
        action="append",
        required=True,
        metavar=("RUN_DIR", "INPUT_DIR"),
    )
    args = ap.parse_args()
    rows = []
    for run, inputs in args.panel:
        for case in csv.DictReader((run / "cases.csv").open()):
            key = {
                k: case[k]
                for k in ("cohort", "target", "query", "receiver", "seed", "arm")
            }
            mode = (
                run
                / "raw"
                / ("{cohort}-{query}-{receiver}-{seed}-{arm}.pdbqt".format(**key))
            )
            row = dict(**key, exit=case["exit"], ligand_top1_rmsd=case.get("top1", ""))
            if case["exit"] == "0" and int(case["reported"]):
                row.update(measure(inputs / key["receiver"], mode))
            rows.append(row)
    with args.output.open("w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=list(dict.fromkeys(k for row in rows for k in row)),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
