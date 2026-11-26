"""Audit official acid-state and flexible-sidechain pools directly from the deposition and PDBQT."""

import argparse
import csv
import itertools
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem

CODES = {"2NNQ": "T4B", "3FR2": "8CA", "1RBP": "RTL", "3FMZ": "2T1"}
ACIDS = {"2NNQ", "3FR2", "3FMZ"}


def xyz(line):
    return np.array([float(line[i : i + 8]) for i in (30, 38, 46)])


def records(text):
    return [row for row in text.splitlines() if row.startswith(("ATOM  ", "HETATM"))]


def protein_atoms(rows):
    return {
        l[12:26]: xyz(l)
        for l in rows
        if l[17:20] != "UNL" and l[12:16].strip()[0] != "H"
    }


def ligand(text):
    # Select residue UNL independently; output flex residue serials reuse ligand numbers.
    labels = next(
        l.split()[2]
        for l in text.splitlines()
        if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX")
    )
    mol = Chem.MolFromSmiles(labels)
    serial = {}
    for l in text.split("BEGIN_RES")[0].splitlines():
        if l.startswith("REMARK SMILES IDX"):
            v = list(map(int, l.split()[3:]))
            serial.update(zip(v[::2], v[1::2], strict=True))
    ordered = {int(l[6:11]): xyz(l) for l in records(text) if l[17:20] == "UNL"}
    positions = np.stack([ordered[serial[i + 1]] for i in range(mol.GetNumAtoms())])
    return mol, positions


def neutral_mol(mol):
    graph = Chem.RWMol(mol)
    for atom in graph.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(False)
    Chem.SanitizeMol(graph)
    return graph.GetMol()


def transform(q, r, graphs):
    x, cx, _ = graphs[q]
    y, cy, _ = graphs[r]
    keys = [
        k
        for k in sorted(x.keys() & y.keys())
        if np.linalg.norm(x[k] - cx) < 14 and np.linalg.norm(y[k] - cy) < 14
    ]
    a, b = np.stack([x[k] for k in keys]), np.stack([y[k] for k in keys])
    u, _s, vt = np.linalg.svd((a - a.mean(0)).T @ (b - b.mean(0)))
    u[:, -1] *= np.sign(np.linalg.det(u @ vt))
    rot = u @ vt
    return lambda p: (p - a.mean(0)) @ rot + b.mean(0)


def matches(crystal, mol, positions, trans, acid):
    reference = trans(crystal.GetConformer().GetPositions())
    neutral = neutral_mol(mol)
    perms = crystal.GetSubstructMatches(
        neutral, uniquify=False, useChirality=True, maxMatches=10001
    )
    assert 0 < len(perms) < 10001
    mapped = np.array([reference[list(m)] for m in perms])

    def rms(arr):
        return float(
            np.linalg.norm((arr - positions).reshape(len(arr), -1), axis=1).min()
            / np.sqrt(len(positions))
        )

    fixed = rms(mapped)
    if not acid:
        return fixed, fixed
    ((_, a, b),) = neutral.GetSubstructMatches(
        Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    )
    swapped = mapped.copy()
    swapped[:, [a, b]] = mapped[:, [b, a]]
    return fixed, min(fixed, rms(swapped))


def write(folder, name, rows):
    with (folder / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def verify(args):
    sys.path.insert(0, str(args.panel.resolve()))
    from verify import reference as derive_crystal

    graphs = {p: derive_crystal(args.fixtures, p, code) for p, code in CODES.items()}
    frames = {
        (q, r): transform(q, r, graphs) for q, r in itertools.product(CODES, repeat=2)
    }
    charges = []
    sidechains = []
    counts = defaultdict(list)
    errors = defaultdict(list)
    population = {}
    target = lambda q: "FABP4" if q in ("2NNQ", "3FR2") else "RBP4"
    for q in sorted(ACIDS):
        raw = {s: (args.acid / f"{q}-{s}.pdbqt").read_text() for s in ("A", "HA")}
        data = {s: ligand(raw[s]) for s in raw}
        neutral, anion = data["HA"][0], data["A"][0]
        assert Chem.MolToSmiles(neutral) == Chem.MolToSmiles(neutral_mol(anion))
        maps = neutral.GetSubstructMatches(
            neutral_mol(anion), uniquify=False, useChirality=True
        )
        err = min(
            np.max(np.linalg.norm(data["HA"][1][list(m)] - data["A"][1], axis=1))
            for m in maps
        )
        assert err == 0
        for state, graph in [("A", anion), ("HA", neutral)]:
            atoms = records(raw[state])
            charge = sum(float(l.split()[-2]) for l in atoms)
            assert abs(charge - (0 if state == "HA" else -1)) < 0.015
            charges.append(
                {
                    "query": q,
                    "state": state,
                    "heavy_atoms": graph.GetNumHeavyAtoms(),
                    "net_formal_charge": Chem.GetFormalCharge(graph),
                    "summed_pdbqt_charge": round(charge, 4),
                    "torsdof": int(
                        next(
                            l.split()[1]
                            for l in raw[state].splitlines()
                            if l.startswith("TORSDOF")
                        )
                    ),
                    "max_input_heavy_shift_A": err,
                }
            )
    for r in CODES:

        def signatures(rows):
            return Counter(
                (l[12:26], tuple(xyz(l)), l.split()[-2], l.split()[-1]) for l in rows
            )

        orig = signatures(records((args.fixtures / r / "receptor.pdbqt").read_text()))
        part = signatures(
            records((args.flex / f"{r}_flex.pdbqt").read_text())
            + records((args.flex / f"{r}_rigid.pdbqt").read_text())
        )
        assert part == orig, r
        population[r] = sum(orig.values())
    for mode, out in [("acid", args.acid), ("flex", args.flex)]:
        schema = (
            ["query", "receiver", "state", "force", "seed", "rank"]
            if mode == "acid"
            else ["query", "receiver", "arm", "seed", "rank"]
        )
        fields = ["rmsd_fixed", "rmsd_resonance"] if mode == "acid" else ["rmsd"]
        recorded = {
            tuple(row[k] for k in schema): [float(row[k]) for k in fields]
            for row in csv.DictReader((out / "modes.csv").open())
        }
        for file in sorted(out.glob("*-*-*-*.pdbqt")):
            parts = file.stem.split("-")
            q, r = parts[:2]
            if mode == "acid":
                state, force, seed = parts[2:]
                arm = state + "_" + force
            else:
                seed, arm = parts[2:]
            trans = frames[q, r]
            values = []
            plain = []
            models = file.read_text().split("ENDMDL")[:-1]
            for rank, m in enumerate(models, 1):
                mol, coords = ligand(m)
                fixed, equiv = matches(graphs[q][2], mol, coords, trans, mode == "acid")
                if mode == "acid":
                    k = (q, r, state, force, seed, str(rank))
                    measured = [fixed, equiv]
                else:
                    k = (q, r, arm, seed, str(rank))
                    measured = [fixed]
                    sides = protein_atoms(records(m))
                    assert bool(sides) == (arm == "flex")
                    if sides:
                        start = protein_atoms(
                            records((out / f"{r}_flex.pdbqt").read_text())
                        )
                        from_query = protein_atoms(
                            records((out / f"{q}_flex.pdbqt").read_text())
                        )
                        assert sides.keys() == start.keys() == from_query.keys()
                        keys = [a for a in sides if a[:4].strip() != "CA"]
                        rms = lambda array: float(
                            np.sqrt(np.mean(np.sum(array**2, axis=1)))
                        )
                        delta = rms(np.stack([sides[a] - start[a] for a in keys]))
                        qdelta = rms(
                            np.stack([sides[a] - trans(from_query[a]) for a in keys])
                        )
                        sidechains.append(
                            {
                                "query": q,
                                "receiver": r,
                                "seed": seed,
                                "rank": rank,
                                "site": keys[0][5:].strip(),
                                "n_sidechain_heavy": len(keys),
                                "rmsd_from_receiver_A": delta,
                                "rmsd_from_transferred_query_A": qdelta,
                            }
                        )
                expected = recorded.pop(k)
                errors[mode].append(
                    max(abs(a - b) for a, b in zip(expected, measured, strict=True))
                )
                values.append(equiv)
                plain.append(fixed)
            for n in (1, 5, 9):
                counts[
                    (mode, target(q), arm, "cross" if q != r else "cognate", n)
                ].append((min(plain[:n]) <= 2, min(values[:n]) <= 2, len(values)))
        assert not recorded
    tables = []
    for (mode, protein, arm, kind, n), cells in sorted(counts.items()):
        tables.append(
            {
                "experiment": mode,
                "target": protein,
                "arm": arm,
                "subset": kind,
                "top_k": n,
                "n_cells": len(cells),
                "hits_fixed": sum(v[0] for v in cells),
                "hits_equiv": sum(v[1] for v in cells),
                "shortest_pool": min(v[2] for v in cells),
            }
        )
    assert max(errors["acid"] + errors["flex"]) < 1e-6
    write(args.out, "charge_checks.csv", charges)
    write(args.out, "threshold_checks.csv", tables)
    write(args.out, "sidechain_rmsd.csv", sidechains)
    summary = {
        "acid_modes": len(errors["acid"]),
        "flex_modes": len(errors["flex"]),
        "acid_max_rmsd_error": max(errors["acid"]),
        "flex_max_rmsd_error": max(errors["flex"]),
        "receptor_atom_records_verified": population,
    }
    (args.out / "audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("panel", "fixtures", "acid", "flex", "out"):
        p.add_argument(arg, type=Path)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    verify(args)
