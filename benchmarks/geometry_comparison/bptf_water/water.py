"""Transfer recurrent BPTF waters from two independent crystal structures."""

import argparse
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "tools"))
sys.path.insert(0, str(HERE.parent))
from extract_independent_receptor_evidence import (
    align_donor,
    aligned_waters,
    structure_chain,
    water_sites,
)
from pose_run import align_with_numbering, load_meta, measure, write_csv
from urllib.request import urlopen

DONORS = ("7K6R", "7M2E")
# Both map by sequence and backbone to the ligand-recognition Asn3007 of 7K6R.
ANCHOR = {"7F5D": 148, "7RWO": 2881}
ARMS = ("dry", "steric", "polar")


def fixed_water(site, serial, residue, polar):
    xyz = np.asarray(site["xyz"])
    atoms = [("O", "OA" if polar else "O", np.zeros(3))]
    if polar:
        # Vina uses bonded HD to type the isotropic oxygen donor/acceptor.
        # Hydrogen orientation here carries no directional water interaction.
        atoms += [
            ("H1", "HD", np.array([0.9572, 0.0, 0.0])),
            ("H2", "HD", np.array([-0.2400, 0.9266, 0.0])),
        ]
    lines = []
    for i, (name, atom_type, offset) in enumerate(atoms):
        x, y, z = xyz + offset
        lines.append(
            f"HETATM{serial + i:5d} {name:4s} HOH W{residue:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    +0.000 {atom_type}"
        )
    return lines


def prepare(inputs, work):
    work.mkdir(parents=True, exist_ok=True)
    for donor in DONORS:
        path = work / (donor + ".pdb")
        if not path.exists():
            path.write_bytes(
                urlopen(
                    f"https://files.rcsb.org/download/{donor}.pdb", timeout=45
                ).read()
            )
    evidence = {}
    for pdb, number in ANCHOR.items():
        _, receptor = structure_chain(inputs / pdb / (pdb + ".pdb"), pdb, "A")
        anchor_key = (number, "")
        anchor = receptor[anchor_key].atoms["ND2"]
        hits, fits = [], {}
        for donor in DONORS:
            path = work / (donor + ".pdb")
            model, protein = structure_chain(path, donor, "A")
            fits[donor] = {"sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            mapping, rot, shift = align_donor(receptor, protein, anchor, fits[donor])
            if mapping[anchor_key] != (3007, "") or protein[(3007, "")].name != "ASN":
                raise ValueError("Unverified pocket residue mapping")
            hits.extend(
                w
                for w in aligned_waters(
                    model, donor, anchor, np.full(3, 16.0), rot, shift
                )
                if w["chain"] == "A" and w["occupancy"] >= 0.8
            )
        dry = (inputs / pdb / "receptor.pdbqt").read_text()
        heavy = np.array(
            [
                [float(l[i : i + 8]) for i in (30, 38, 46)]
                for l in dry.splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
                and l[77:].strip() not in ("H", "HD")
            ]
        )
        selected = []
        for s in water_sites(hits):
            xyz = np.array(s["xyz"])
            if (
                s["independent_pdb_count"] == 2
                and np.linalg.norm(xyz - anchor) <= 8
                and np.linalg.norm(heavy - xyz, axis=1).min() >= 2.5
            ):
                selected.append(s)
        if not selected:
            raise ValueError("No independently recurrent water in pocket")
        for i, s in enumerate(selected):
            if any(
                np.linalg.norm(np.subtract(s["xyz"], t["xyz"])) < 2.4
                for t in selected[:i]
            ):
                raise ValueError("Incompatible adjacent water sites")
        evidence[pdb] = dict(fits=fits, anchor=anchor.tolist(), sites=selected)
        for arm in ARMS:
            out = work / arm / pdb
            out.mkdir(parents=True, exist_ok=True)
            ligand = out / "ligand.pdbqt"
            ligand.write_bytes((inputs / pdb / "ligand.pdbqt").read_bytes())
            extra = []
            if arm != "dry":
                for i, s in enumerate(selected):
                    extra.extend(fixed_water(s, 9001 + 3 * i, 800 + i, arm == "polar"))
            (out / "receptor.pdbqt").write_text(
                dry.rstrip() + "\n" + "".join(l + "\n" for l in extra)
            )
    (work / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    a = p.parse_args()
    prepare(a.inputs, a.work)
    meta = {pdb: load_meta(a.inputs, pdb) for pdb in ANCHOR}
    rows, modes = [], []
    for cap, seeds in ((48000, (84551, 84559, 84571)), (120000, (93701, 93719, 93731))):
        for q in ANCHOR:
            for r in ANCHOR:
                rot, tr, rms, n, offset = np.eye(3), np.zeros(3), 0.0, 0, 0
                if q != r:
                    rot, tr, rms, n, offset = align_with_numbering(meta[q], meta[r])
                ref = np.asarray(meta[q]["mapped"]) @ rot + tr
                for i, seed in enumerate(seeds):
                    for arm in ARMS[i:] + ARMS[:i]:
                        out = a.work / arm / "poses"
                        out.mkdir(exist_ok=True)
                        row, poses = measure(
                            (
                                a.work / arm,
                                out,
                                a.binary.resolve(),
                                meta[q],
                                meta[r],
                                ref,
                                rms,
                                n,
                                offset,
                                seed,
                                arm,
                                cap,
                            )
                        )
                        rows.append(row)
                        modes.extend(poses)
    write_csv(a.work / "cells.csv", rows)
    write_csv(a.work / "modes.csv", modes)


if __name__ == "__main__":
    main()
