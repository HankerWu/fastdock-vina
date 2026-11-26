"""Compute fixed-pool scores, matched chemical contrasts, crystal transfer and Ki direction."""

import argparse
import itertools
import math
from pathlib import Path
from statistics import median

import numpy as np
from rdkit import Chem

from experiment import ARMS, ROOT, SEEDS, aspsite, table, write
from make_ion_pair_anchors import ligand_sites, site_points, xyz
from trypsin import RECEPTORS, SERIES

RT = 0.00198720426 * 298.15
ODDS = math.expm1(2.0)
METHODS = ("release", "pair_rerank", "competition")


def point_error(coords, targets):
    # Equivalent amidinium nitrogens and carboxylate oxygens give two assignments.
    return min(float(np.square(coords - t).sum() / 2) for t in (targets, targets[::-1]))


def site_measure(group, targets, coords):
    return point_error(np.array([coords[i] for i in group[1:]]), targets)


def pose_models(path):
    blocks = [b.splitlines() for b in path.read_text().split("ENDMDL") if b.strip()]
    for lines in blocks:
        affinity = float(
            next(l.split()[3] for l in lines if l.startswith("REMARK VINA RESULT:"))
        )
        coords = {
            int(l[6:11]): np.array(xyz(l))
            for l in lines
            if l.startswith(("ATOM", "HETATM"))
        }
        yield affinity, coords


def ring_maps(ligand):
    lines = ligand.read_text().splitlines()
    from score_physical_pool import serial_index_pairs

    serials = {i - 1: s for i, s in serial_index_pairs(lines, "REMARK SMILES IDX ")}
    smiles = next(
        l.split()[2]
        for l in lines
        if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX")
    )
    mol = Chem.MolFromSmiles(smiles)
    # A carbocation resonance representation changes primary N valence. Allow either N bond order.
    p = Chem.MolFromSmarts("[c]1(-[C](~[#7])~[#7])[c][c][c][c][c]1")
    maps = mol.GetSubstructMatches(p, uniquify=False)
    if len(maps) != 4:
        raise ValueError(f"Expected four benzamidinium graph mappings, got {len(maps)}")
    return [[serials[i] for i in m] for m in maps]


def ben_crystal(path):
    atoms = {
        l[12:16].strip(): np.array(xyz(l))
        for l in path.read_text().splitlines()
        if l.startswith("HETATM")
        and l[17:20] == "BEN"
        and l[21] == "A"
        and l[16] in " A"
    }
    central = atoms["C"]
    rings = ["C" + str(i) for i in range(1, 7)]
    start = int(np.argmin([np.linalg.norm(central - atoms[x]) for x in rings]))
    seq = rings[start:] + rings[:start]
    return np.array([atoms[n] for n in [seq[0], "C", "N1", "N2", *seq[1:]]])


def external_sites(raw, fjs):
    out = []
    for name, comp, path in (
        ("3PTB", "BEN", raw / "3PTB.pdb"),
        ("1BTY", "BEN", raw / "1BTY.pdb"),
        ("5MO0", "BEN", raw / "5MO0.pdb"),
        ("1FAX", "DX9", raw / "1FAX.pdb"),
        ("1FJS", "Z34", fjs / "1FJS.pdb"),
    ):
        ca = aspsite(path, "A")
        targets = site_points(ca[0], ca[1:])
        records = [
            l
            for l in path.read_text().splitlines()
            if l.startswith("HETATM")
            and l[17:20] == comp
            and l[21] == "A"
            and l[16] in " A"
        ]
        cs = {
            l[12:16].strip(): np.array(xyz(l))
            for l in records
            if l[76:78].strip() == "C"
        }
        ns = {
            l[12:16].strip(): np.array(xyz(l))
            for l in records
            if l[76:78].strip() == "N"
        }
        groups = [
            (c, [(n, v) for n, v in ns.items() if np.linalg.norm(v - vc) < 1.60])
            for c, vc in cs.items()
        ]
        matches = [
            (point_error(np.array([v for n, v in pair]), targets), c, pair)
            for c, pair in groups
            if len(pair) == 2
        ]
        err, central, pair = min(matches, key=lambda m: m[0])
        a = np.array([v for _, v in pair])
        b = ca[1:]
        distances = min(
            [sorted(np.linalg.norm(a - o, axis=1)) for o in (b, b[::-1])], key=max
        )
        out.append(
            dict(
                pdb=name,
                ligand=comp,
                C=central,
                names=",".join(n for n, v in pair),
                site_rms=math.sqrt(err),
                short_NO=distances[0],
                long_NO=distances[1],
            )
        )
    return out


def panel(target, ligands, work, raw, orig, seeds=SEEDS, arms=ARMS):
    modes, cells = [], []
    for rec, lig, seed, arm in itertools.product(orig, ligands, seeds, arms):
        inputs = orig[rec]
        if target == "fxa":
            input_lig = inputs / f"{lig}.pdbqt"
            structure = inputs / f"{rec}.pdb"
            path = (
                inputs / f"{lig}_{seed}.pdbqt"
                if arm == "native" and seed in SEEDS
                else work / f"{rec}_{lig}_{seed}_{arm}.pdbqt"
            )
        else:
            input_lig = inputs / f"{lig}.pdbqt"
            structure = raw / f"{rec}.pdb"
            path = work / f"{rec}_{lig}_{seed}_{arm}.pdbqt"
        groups, acids = ligand_sites(input_lig)
        ca = aspsite(structure, "A")
        targets = site_points(ca[0], ca[1:])
        maps = ring_maps(input_lig) if target == "trypsin" and lig == "H" else None
        crystal = ben_crystal(structure) if maps else None
        pool = []
        for model, (score, coords) in enumerate(pose_models(path), 1):
            errors = [site_measure(g, targets, coords) for g in groups]
            bound = min(errors)
            match = math.exp(
                -0.5 * bound
            )  # squared site RMS displacement in Å², width 1 Å
            intra_err = []
            for acid in acids:
                s = np.array([coords[i] for i in acid])
                self_target = site_points(s[0], s[1:])
                intra_err.extend(site_measure(g, self_target, coords) for g in groups)
            self_contact = math.exp(-0.5 * min(intra_err)) if intra_err else 0.0
            selected = groups[int(np.argmin(errors))]
            n = np.asarray([coords[i] for i in selected[1:]])
            contact = min(
                max(np.linalg.norm(n - o, axis=1)) for o in (ca[1:], ca[1:][::-1])
            )
            pose_rms = (
                min(
                    np.sqrt(
                        np.sum(np.square(np.array([coords[i] for i in m]) - crystal))
                        / len(m)
                    )
                    for m in maps
                )
                if maps
                else ""
            )
            group = dict(target=target, receptor=rec, ligand=lig, seed=seed, arm=arm)
            pool.append(
                dict(
                    **group,
                    model=model,
                    release=score,
                    pair_rerank=score - RT * math.log1p(ODDS * match),
                    competition=score
                    - RT * math.log1p(ODDS * match / (1 + ODDS * self_contact)),
                    bound_rms=math.sqrt(bound),
                    self_contact=self_contact,
                    largest_NO=contact,
                    heavy_rms=pose_rms,
                )
            )
        modes.extend(pool)
        for method in METHODS:
            top = min(pool, key=lambda row: row[method])
            cells.append(
                dict(
                    **group,
                    method=method,
                    model=top["model"],
                    score=top[method],
                    bound_rms=top["bound_rms"],
                    bidentate=int(top["largest_NO"] <= 3.5),
                    oracle_bidentate=int(any(r["largest_NO"] <= 3.5 for r in pool)),
                    heavy_rms=top["heavy_rms"],
                    oracle_rms=min(r["heavy_rms"] for r in pool),
                )
            )
    return modes, cells


def compare(cells, labels):
    pairs = []
    for (t, arm, method), rows in itertools.groupby(
        sorted(cells, key=lambda r: (r["target"], r["arm"], r["method"])),
        key=lambda r: (r["target"], r["arm"], r["method"]),
    ):
        by = {(r["receptor"], r["seed"], r["ligand"]): r for r in rows}
        panels = sorted(set((r, s) for r, s, l in by))
        for a, b in itertools.combinations(sorted(labels[t], key=labels[t].get), 2):
            ka, kb = labels[t][a], labels[t][b]
            if ka == kb:
                continue
            gaps = [by[r, s, b]["score"] - by[r, s, a]["score"] for r, s in panels]
            gap = round(median(gaps), 7)
            pairs.append(
                dict(
                    target=t,
                    arm=arm,
                    method=method,
                    strong=a,
                    weak=b,
                    ratio=kb / ka,
                    experimental_gap=RT * math.log(kb / ka),
                    predicted_gap=gap,
                    credit=1 if gap > 0 else 0 if gap < 0 else 0.5,
                    min_gap=min(gaps),
                    max_gap=max(gaps),
                )
            )
    return pairs


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("work", "fjs", "p16", "out"):
        p.add_argument(key, type=Path)
    p.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    a = p.parse_args()
    a.out.mkdir(exist_ok=True, parents=True)
    metadata = [
        r
        for r in table(ROOT / "benchmarks/affinity_partition/fxa_series.csv")
        if int(r["included"])
    ]
    labels = {
        "fxa": {r["ligand"]: float(r["ki_nm"]) for r in metadata},
        "trypsin": {k: 100 / v for k, (s, v) in SERIES.items()},
    }
    sizes = {r["ligand"]: (r["heavy"], r["charge"]) for r in metadata}
    all_modes, all_cells = [], []
    for target, inputs in (
        ("fxa", {"1FJS": a.fjs, "2P16": a.p16}),
        ("trypsin", {k: a.work / "trypsin" for k in RECEPTORS}),
    ):
        m, c = panel(
            target, labels[target], a.work / target, a.work / "raw", inputs, a.seeds
        )
        all_modes += m
        all_cells += c
    pairs = compare(all_cells, labels)
    write(a.out / "modes.csv", all_modes)
    write(a.out / "cells.csv", all_cells)
    write(a.out / "pairs.csv", pairs)
    write(a.out / "crystals.csv", external_sites(a.work / "raw", a.fjs))
    summary = []
    for t, arm, method in itertools.product(labels, ARMS, METHODS):
        group = [
            r
            for r in all_cells
            if (r["target"], r["arm"], r["method"]) == (t, arm, method)
        ]
        edges = [
            r
            for r in pairs
            if (r["target"], r["arm"], r["method"]) == (t, arm, method)
            and r["ratio"] >= 3
        ]
        matched = [
            r for r in edges if t == "fxa" and sizes[r["strong"]] == sizes[r["weak"]]
        ]
        summary.append(
            dict(
                target=t,
                arm=arm,
                method=method,
                cells=len(group),
                bidentate=sum(r["bidentate"] for r in group),
                oracle_bidentate=sum(r["oracle_bidentate"] for r in group),
                unique=len(edges),
                credit=sum(r["credit"] for r in edges),
                equal_size_charge=len(matched),
                matched_credit=sum(r["credit"] for r in matched),
            )
        )
    write(a.out / "summary.csv", summary)
