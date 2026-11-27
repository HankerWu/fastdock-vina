"""Validate two-state algebra and replay one terminal-identity change in the native scorer."""

import argparse
import csv
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from model import RT, AmideReceptor, atoms, marginal
from scipy.spatial.transform import Rotation

CASES = [
    ("3TUH", 3000, 43103, 2, "ASN A 106 "),
    ("6I62", 3000, 43117, 8, "GLN A 426 "),
    ("6I62", 3000, 43103, 3, "ASN A 346 "),
    ("2XJG", 12000, 43103, 1, "ASN A  51 "),
]


def check_partition():
    for a in (-500.0, -2.0, 0.0, 3.0, 500.0):
        for dl in (-3.0, 0.0, 4.0):
            result = marginal(a, dl)
            # Conditioning the reference identity on the other state changes the baseline by dl.
            assert abs(result - dl - marginal(-a, -dl)) < 2.0e-12
            assert min(0.0, dl) - 1.0e-12 <= result <= max(0.0, dl) + 1.0e-12
            assert (
                abs(
                    result
                    + RT
                    * np.log(
                        (1 - 1 / (1 + np.exp(np.clip(a / RT, -700, 700))))
                        + np.exp(-dl / RT) / (1 + np.exp(np.clip(a / RT, -700, 700)))
                    )
                )
                < 2.0e-12
            )


def flip(pdbqt, residue):
    group = {i: a for i, a in atoms(pdbqt).items() if a.residue == residue}
    names = {a.name: i for i, a in group.items()}
    o, n, c = ("OD1", "ND2", "CG") if residue[:3] == "ASN" else ("OE1", "NE2", "CD")
    oi, ni, ci = [names[k] for k in (o, n, c)]
    ox, nx, cx = [group[i].xyz for i in (oi, ni, ci)]
    a, b = nx - cx, ox - cx
    axis = np.cross(a, b)
    angle = np.arctan2(np.linalg.norm(axis), np.dot(a, b))
    rot = Rotation.from_rotvec(angle * axis / np.linalg.norm(axis))
    new = {oi: nx, ni: ox}
    # Move attached H with the newly located amide nitrogen. The native heavy-atom
    # Vina potential uses these polar hydrogens for donor identity only.
    hs = [
        i
        for i, at in group.items()
        if at.kind == "HD" and np.linalg.norm(at.xyz - nx) < 1.25
    ]
    assert len(hs) == 2
    for h in hs:
        new[h] = ox + rot.apply(group[h].xyz - nx)
    out = []
    for line in pdbqt.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and int(line[6:11]) in new:
            line = (
                line[:30]
                + "".join(f"{v:8.3f}" for v in new[int(line[6:11])])
                + line[54:]
            )
        out.append(line)
    return "\n".join(out) + "\n"


def intermolecular(binary, receptor, ligand, pose):
    cmd = [
        str(binary),
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--pose_coordinates",
        str(pose),
        "--score_only",
        "--net_polar",
        "combined",
        "--cpu",
        "1",
        "--seed",
        "42",
        "--verbosity",
        "0",
    ]
    text = subprocess.check_output(cmd, text=True)
    line = next(
        line for line in text.splitlines() if line.startswith("Net polar native:")
    )
    return float(line.split(":", 1)[1].split()[0])


def replay(work, binary, output):
    check_partition()
    rows = []
    for q, cap, seed, rank, residue in CASES:
        folder = work / q
        receptor = folder / "receptor.pdbqt"
        ligand = folder / "ligand.pdbqt"
        before = receptor.read_text()
        block = (
            (work / "runs" / f"{q}_{q}_{cap}_{seed}_pro.pdbqt")
            .read_text()
            .split("ENDMDL")[rank - 1]
        )
        terms = AmideReceptor(before).score(ligand.read_text(), block, detail=True)
        _, _, delta, _ = next(s for s in terms["sites"] if s[0] == residue)
        with TemporaryDirectory() as temp:
            path = Path(temp)
            alt = path / "flip.pdbqt"
            pose = path / "pose.pdbqt"
            alt.write_text(flip(before, residue))
            pose.write_text(
                "\n".join(l for l in block.splitlines() if not l.startswith("MODEL"))
                + "\n"
            )
            e0 = intermolecular(binary, receptor, ligand, pose)
            e1 = intermolecular(binary, alt, ligand, pose)
        rows.append(
            {
                "query": q,
                "cap": cap,
                "seed": seed,
                "model": rank,
                "residue": residue,
                "analytic_delta": delta,
                "native_deposited": e0,
                "native_flipped": e1,
                "native_delta": e1 - e0,
                "residual": delta - (e1 - e0),
            }
        )
    with output.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    # Table interpolation in the production pair potential and printed coordinates/scores can differ slightly.
    assert max(abs(r["residual"]) for r in rows) < 0.01, rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    replay(args.work, args.binary.resolve(), args.output)
