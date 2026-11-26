"""Audit output energies, contacts, independent pose RMSD and distinct assay pairs."""

import argparse
import csv
import itertools
import math
from collections import defaultdict
from pathlib import Path
from statistics import median
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def read_coords(lines):
    coords = {}
    elements = {}
    for l in lines:
        if l[:6].strip() in ("ATOM", "HETATM"):
            i = int(l[6:11])
            coords[i] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            elements[i] = l.split()[-1]
    return coords, elements


def amidine(lig):
    coords, e = read_coords(lig.read_text().splitlines())
    for c, t in e.items():
        if t != "C":
            continue
        ns = [
            n
            for n, t in e.items()
            if t == "N" and np.linalg.norm(coords[n] - coords[c]) < 1.65
        ]
        if len(ns) != 2:
            continue
        hyd = [
            sum(
                t == "HD" and np.linalg.norm(coords[h] - coords[n]) < 1.24
                for h, t in e.items()
            )
            for n in ns
        ]
        if hyd == [2, 2] and any(
            t == "A" and np.linalg.norm(coords[a] - coords[c]) < 1.65
            for a, t in e.items()
        ):
            return ns
    raise ValueError(
        f"Physical protonation/topology cannot identify one primary aryl amidine in {lig}"
    )


def match_benzamidine(prep, crys):
    a, types = read_coords(prep.read_text().splitlines())
    source = {i: ("N" if t == "N" else "C") for i, t in types.items() if t != "HD"}
    lines = [
        l
        for l in crys
        if l.startswith("HETATM")
        and l[17:20] == "BEN"
        and l[21] == "A"
        and l[16] in " A"
        and l[76:78].strip() in ("C", "N")
    ]
    b = {
        int(l[6:11]): np.array([float(l[k : k + 8]) for k in (30, 38, 46)])
        for l in lines
    }
    target = {int(l[6:11]): l[76:78].strip() for l in lines}

    def edges(nodes, coords):
        return {
            frozenset((i, j))
            for i, j in itertools.combinations(nodes, 2)
            if np.linalg.norm(coords[i] - coords[j]) < 1.7
        }

    src_edges = edges(source, a)
    dst_edges = edges(target, b)

    def group_by_element_and_degree(nodes, bonds):
        groups = defaultdict(list)
        for i, symbol in nodes.items():
            groups[symbol, sum(i in edge for edge in bonds)].append(i)
        return groups

    groups_a = group_by_element_and_degree(source, src_edges)
    groups_b = group_by_element_and_degree(target, dst_edges)
    ids = [j for group in groups_a.values() for j in group]
    mats = []
    for assignment in itertools.product(
        *(itertools.permutations(groups_b[k]) for k in groups_a)
    ):
        mapping = dict(zip(ids, itertools.chain.from_iterable(assignment), strict=True))
        if {frozenset(mapping[i] for i in e) for e in src_edges} == dst_edges:
            mats.append(mapping)
    if len(mats) != 4:
        raise ValueError("Benzamidinium automorphism mismatch")
    return mats, b


def verify(work, fjs, p16, results):
    groups = defaultdict(list)
    observed = defaultdict(dict)
    for r in rows(results / "modes.csv"):
        groups[
            r["target"], r["receptor"], r["ligand"], int(r["seed"]), r["arm"]
        ].append(r)
    reference = {}
    ns_cache = {}
    mats_cache = {}
    nposes = 0
    maxd = 0
    maxe = 0
    maxr = 0
    for (target, rec, lig, seed, arm), data in groups.items():
        wd = fjs if rec == "1FJS" else p16 if rec == "2P16" else work / "trypsin"
        mol = wd / f"{lig}.pdbqt"
        if (rec, lig) not in ns_cache:
            ns_cache[rec, lig] = amidine(mol)
        if rec not in reference:
            lines = (
                ((work / "raw" if target == "trypsin" else wd) / f"{rec}.pdb")
                .read_text()
                .splitlines()
            )
            os = [
                l
                for l in lines
                if l.startswith("ATOM")
                and l[12:16].strip() in ("OD1", "OD2")
                and l[22:26].strip() == "189"
                and l[17:20] == "ASP"
                and l[21] == "A"
                and l[16] in " A"
            ]
            os = sorted(os, key=lambda l: l[12:16])
            assert len(os) == 2
            reference[rec] = np.array(
                [[float(l[k : k + 8]) for k in (30, 38, 46)] for l in os]
            )
        if lig == "H" and rec not in mats_cache:
            mats_cache[rec] = match_benzamidine(
                mol, (work / "raw" / f"{rec}.pdb").read_text().splitlines()
            )
        if target == "fxa" and arm == "native" and seed in (271828, 314159):
            path = wd / f"{lig}_{seed}.pdbqt"
        else:
            path = work / target / f"{rec}_{lig}_{seed}_{arm}.pdbqt"
        blocks = path.read_text().split("MODEL")[1:]
        assert len(blocks) == len(data)
        for block, r in zip(blocks, data, strict=True):
            lines = block.splitlines()
            x, _ = read_coords(lines)
            score = float(
                next(
                    l[len("REMARK VINA RESULT:") :].split()[0]
                    for l in lines
                    if l.startswith("REMARK VINA RESULT:")
                )
            )
            if int(r["model"]) == 1:
                observed[target, arm, lig][rec, seed] = score
            maxe = max(maxe, abs(score - float(r["release"])))
            ns = [x[k] for k in ns_cache[rec, lig]]
            ds = [[math.dist(n, o) for o in reference[rec]] for n in ns]
            best = min(max(ds[0][0], ds[1][1]), max(ds[0][1], ds[1][0]))
            maxd = max(maxd, abs(best - float(r["largest_NO"])))
            if lig == "H":
                maps, template = mats_cache[rec]
                rms = min(
                    math.sqrt(
                        sum(math.dist(x[i], template[j]) ** 2 for i, j in m.items())
                        / len(m)
                    )
                    for m in maps
                )
                maxr = max(maxr, abs(rms - float(r["heavy_rms"])))
            nposes += 1
    assert maxr < 1e-10 and maxe < 1e-12 and maxd < 1e-10
    checked = 0
    for row in rows(results / "pairs.csv"):
        if row["method"] != "release":
            continue
        d = lambda lig: observed[row["target"], row["arm"], lig]
        a, b = d(row["strong"]), d(row["weak"])
        delta = median(b[k] - v for k, v in a.items())
        assert abs(delta - float(row["predicted_gap"])) < 1e-6
        checked += 1
    return f"{len(groups)} raw pose files; {nposes} scores/NO geometries; {checked} unique-pair rows reconstructed. Max score/NO/RMSD discrepancies {maxe:.3g}/{maxd:.3g}/{maxr:.3g}."


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "fjs", "p16", "results"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    print(verify(a.work, a.fjs, a.p16, a.results))
