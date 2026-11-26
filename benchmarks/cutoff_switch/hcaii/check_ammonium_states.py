"""Independently check heavy-atom sign and receptor-frame RMSD in the two-state runs."""

import argparse
import csv
import itertools
import json
import re
from collections import defaultdict
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign


def read_points(text):
    return {
        int(line[6:11]): np.array([float(line[j : j + 8]) for j in (30, 38, 46)])
        for line in text.splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    }


def ccd_vectors(work, lig, meta):
    code = meta["ccd"]
    ccd = gemmi.cif.read_file(str(work / "raw" / f"{code}.cif")).sole_block()
    elem = dict(ccd.find("_chem_comp_atom.", ["atom_id", "type_symbol"]))
    adj = defaultdict(list)
    for a, b, order in ccd.find(
        "_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]
    ):
        if "H" not in (elem[a], elem[b]):
            adj[a].append((b, order))
            adj[b].append((a, order))
    # Third-party raw CCD topology: tertiary N with three saturated carbon neighbours.
    (c,) = [
        name
        for name in elem
        if elem[name] == "N"
        and len(adj[name]) == 3
        and all(
            elem[n] == "C" and all(o == "SING" for _, o in adj[n]) for n, _ in adj[name]
        )
    ]
    neighbours = sorted(n for n, _ in adj[c])
    names = [c, *neighbours]
    serials = [meta["qt_serial"][meta["names"].index(n)] for n in names]
    lines = (work / "raw" / f"{lig}.pdb").read_text().splitlines()
    coords = {
        line[12:16].strip(): np.array([float(line[j : j + 8]) for j in (30, 38, 46)])
        for line in lines
        if line.startswith("HETATM")
        and line[17:20] == code
        and line[21] == "A"
        and int(line[22:26]) == 303
        and line[16] in " A"
    }
    return serials, handed([coords[n] for n in names])


def handed(points):
    center, a, b, c = points
    return np.sign(np.dot(a - center, np.cross(b - center, c - center)))


def verify(work, folder):
    meta = json.loads((work / "prepared" / "metadata.json").read_text())
    labels = {
        (r["receptor"], r["ligand"], r["seed"], r["slot"], int(r["model"])): r
        for r in csv.DictReader((folder / "candidates.csv").open())
    }
    observed = {}
    peak_rms_error = 0.0
    peak_distance_change = 0.0
    for lig, data in meta.items():
        serials, reference_sign = ccd_vectors(work, lig, data)
        original = read_points((work / "prepared" / f"{lig}_ligand.pdbqt").read_text())
        mirror = read_points((folder / "inputs" / f"{lig}_reflected.pdbqt").read_text())
        assert original.keys() == mirror.keys()
        delta = max(
            abs(
                np.linalg.norm(original[a] - original[b])
                - np.linalg.norm(mirror[a] - mirror[b])
            )
            for a, b in itertools.combinations(original, 2)
        )
        peak_distance_change = max(delta, peak_distance_change)
        assert handed([original[s] for s in serials]) == -handed(
            [mirror[s] for s in serials]
        )
        graph = Chem.SDMolSupplier(str(work / "prepared" / f"{lig}_graph.sdf"))[0]
        for path in sorted(folder.glob(f"*_{lig}_*.pdbqt")):
            rec, _, seed, slot = path.stem.split("_", maxsplit=3)
            reference = Chem.Mol(graph)
            for idx, pt in enumerate(data["references"][rec]["coords"]):
                reference.GetConformer().SetAtomPosition(idx, pt)
            for number, text in re.findall(
                r"^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL$",
                path.read_text(),
                re.MULTILINE | re.DOTALL,
            ):
                key = (rec, lig, seed, slot, int(number))
                row = labels[key]
                points = read_points(text)
                sign = handed([points[s] for s in serials])
                assert sign == handed(
                    [
                        (mirror if slot == "second_reflected" else original)[s]
                        for s in serials
                    ]
                )
                same = sign == reference_sign
                assert same == bool(int(row["same_sign_as_bound"]))
                probe = Chem.Mol(graph)
                for idx, s in enumerate(data["qt_serial"]):
                    probe.GetConformer().SetAtomPosition(idx, points[s])
                rms = rdMolAlign.CalcRMS(
                    probe,
                    reference,
                    maxMatches=1000000,
                    symmetrizeConjugatedTerminalGroups=True,
                )
                peak_rms_error = max(abs(rms - float(row["rmsd"])), peak_rms_error)
                observed[key] = (rms, same)
    assert (
        set(observed) == set(labels)
        and peak_rms_error < 1e-6
        and peak_distance_change < 1e-10
    )
    totals = defaultdict(lambda: [0, 0, 0, 0])
    for row in csv.DictReader((folder / "cells.csv").open()):
        rec, lig, seed, arm = (row[n] for n in ("receptor", "ligand", "seed", "arm"))
        second = "second_original" if arm == "original_only" else "second_reflected"
        pool = {
            k: val
            for k, val in observed.items()
            if k[:3] == (rec, lig, seed) and k[3] in ("common", second)
        }
        key = (rec, lig, seed, row["selected_slot"], int(row["selected_model"]))
        assert key in pool
        r, same = pool[key]
        counts = [
            r <= 2,
            any(v[0] <= 2 for v in pool.values()),
            r <= 2 and same,
            any(v[0] <= 2 and v[1] for v in pool.values()),
        ]
        totals[arm] = [i + int(j) for i, j in zip(totals[arm], counts, strict=True)]
    print(
        f"poses={len(labels)} maximum_RMSD_delta={peak_rms_error:.9g} max_intramolecular_distance_change={peak_distance_change:.9g}"
    )
    for arm, counts in sorted(totals.items()):
        print(arm, "RMSD top/oracle and sign-compatible top/oracle", counts)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("runs", type=Path)
    args = parser.parse_args()
    verify(args.work, args.runs)
