"""Audit locked ligand trees, deposited CA crystal poses and equal-budget rankings."""

import argparse
import csv
import itertools
import math
import re
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmarks/anchor_search"))
from verify_hcaii import audited_reference, backbone


def atom_line(line):
    return line.startswith(("ATOM  ", "HETATM"))


def atoms(lines):
    return {
        int(l[6:11]): np.array([float(l[k : k + 8]) for k in (30, 38, 46)])
        for l in lines
        if atom_line(l)
    }


def metadata(lines, tag):
    result = {}
    for l in lines:
        if l.startswith(tag):
            seq = list(map(int, l[len(tag) :].split()))
            result.update(zip(seq[::2], seq[1::2], strict=True))
    return result


def geometry(lines):
    smiles = next(
        l.split()[2]
        for l in lines
        if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX ")
    )
    mol = Chem.MolFromSmiles(smiles)
    index = metadata(lines, "REMARK SMILES IDX ")
    serial = [index[i + 1] for i in range(mol.GetNumAtoms())]
    edges = {
        tuple(sorted((serial[b.GetBeginAtomIdx()], serial[b.GetEndAtomIdx()])))
        for b in mol.GetBonds()
    }
    for l in lines:
        if l.startswith("REMARK H PARENT "):
            nums = list(map(int, l.split()[3:]))
            edges.update(
                tuple(sorted((index[i], j)))
                for i, j in zip(nums[::2], nums[1::2], strict=True)
            )
    ranks = Chem.CanonicalRankAtoms(mol, breakTies=False)
    sites = []
    for o, c, n in mol.GetSubstructMatches(Chem.MolFromSmarts("O=[C;X3]-[N;X3;!+]")):
        cn = mol.GetBondBetweenAtoms(c, n)
        alkyl = [
            x.GetIdx() for x in mol.GetAtomWithIdx(n).GetNeighbors() if x.GetIdx() != c
        ]
        if (
            cn.IsInRing()
            or not alkyl
            or any(mol.GetAtomWithIdx(x).GetAtomicNum() != 6 for x in alkyl)
        ):
            continue
        if len({ranks[x] for x in alkyl}) == len(alkyl):
            sites.append(
                tuple(serial[x] for x in (o, c, n, max(alkyl, key=ranks.__getitem__)))
            )
    assert len(sites) <= 1
    return mol, serial, edges, sites[0] if sites else None


def angle(coords, ids):
    a, b, c, d = (coords[i] for i in ids)
    axis = c - b
    n1 = np.cross(b - a, axis)
    n2 = np.cross(axis, d - c)
    return math.degrees(
        math.atan2(
            np.dot(np.cross(n1, n2), axis) / np.linalg.norm(axis), np.dot(n1, n2)
        )
    )


def angular_error(a, b):
    return abs((a - b + 180) % 360 - 180)


def torsion_tree(lines):
    stack = [None]
    parent = {}
    owners = {}
    order = []
    opened = False
    for l in lines:
        t = l.split()
        if t[0] == "ROOT":
            opened = True
        if t[0] == "ENDROOT":
            opened = False
        if t[0] == "BRANCH":
            edge = tuple(map(int, t[1:]))
            assert owners[edge[0]] == stack[-1]
            parent[edge] = stack[-1]
            stack.append(edge)
            order.append(edge)
        if t[0] == "ENDBRANCH":
            assert stack.pop() == tuple(map(int, t[1:]))
        if atom_line(l):
            i = int(l[6:11])
            assert (opened or stack[-1]) and i not in owners
            owners[i] = stack[-1]
    assert stack == [None] and not opened
    for edge in order:
        assert owners[edge[1]] == edge
    (count,) = [int(l.split()[1]) for l in lines if l.startswith("TORSDOF ")]
    return parent, owners, order, count


def audit_inputs(base, matrix):
    maximum_geometry = maximum_dihedral = maximum_target = 0.0
    nstates = already_fixed = 0
    pyramidal = []
    for path in sorted((matrix / "inputs").glob("*_free.pdbqt")):
        panel, lig, state, _ = path.stem.split("_")
        if panel == "ca" and state == "mirror":
            original = (
                base / "fdv-map-states-declared/inputs" / f"{lig}_reflected.pdbqt"
            )
        elif panel == "fxa":
            original = (
                base
                / ("fdv-ionpair/repair" if lig == "CHEMBL369480" else "fdv-fxa1")
                / f"{lig}.pdbqt"
            )
        else:
            original = (
                base
                / ("fdv-second" if panel == "ca" else "fdv-brd4")
                / "prepared"
                / f"{lig}_ligand.pdbqt"
            )
        assert path.read_bytes() == original.read_bytes(), path
        native = path.read_text().splitlines()
        mol, serial, edges, site = geometry(native)
        trusted = {
            tuple(sorted(map(int, l.split())))
            for l in path.with_name(path.stem[:-5] + ".bonds").read_text().splitlines()
        }
        assert edges == trusted
        p0, owners0, order0, t0 = torsion_tree(native)
        if site is None:
            assert not path.with_name(path.name.replace("_free.", "_near.")).exists()
            continue
        o, c, n, tip = site
        target = next((pair for pair in order0 if set(pair) == {c, n}), None)
        already_fixed += target is None
        start = atoms(native)
        adjacent = defaultdict(set)
        for a, b in edges:
            adjacent[a].add(b)
            adjacent[b].add(a)
        side = {n}
        for _ in range(len(start)):
            side.update(set.union(*(adjacent[i] for i in side)) - {c})
        assert c not in side and side | (set(start) - side) == set(start)
        phi0 = angle(start, site)
        for form in ("near", "far"):
            result = (
                path.with_name(path.name.replace("_free.", "_" + form + "."))
                .read_text()
                .splitlines()
            )
            mol1, serial1, edges1, site1 = geometry(result)
            assert (
                Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol)
                and serial1 == serial
                and edges1 == edges
                and site1 == site
            )
            assert sorted(l[:30] + l[54:] for l in native if atom_line(l)) == sorted(
                l[:30] + l[54:] for l in result if atom_line(l)
            )
            p, owners, order, t = torsion_tree(result)
            assert order == [e for e in order0 if e != target] and t == t0 - int(
                target is not None
            )

            def contracted(v):
                return p0[target] if target is not None and v == target else v

            assert all(owners[i] == contracted(v) for i, v in owners0.items())
            assert all(p[e] == contracted(p0[e]) for e in order)
            xyz = atoms(result)
            assert all(
                np.array_equal(xyz[i], start[i]) for i in (set(start) - side) | {n}
            )
            for group in (side | {c}, (set(start) - side) | {n}):
                for a, b in itertools.combinations(group, 2):
                    delta = abs(
                        np.linalg.norm(xyz[a] - xyz[b])
                        - np.linalg.norm(start[a] - start[b])
                    )
                    maximum_geometry = max(maximum_geometry, delta)
                    assert delta < 0.00175
            for b, c2 in edges:
                if {b, c2} == {site[1], site[2]}:
                    continue
                for a, d in itertools.product(adjacent[b] - {c2}, adjacent[c2] - {b}):
                    if a == d:
                        continue
                    # A collinear triple has no well-defined torsion.
                    if (
                        min(
                            np.linalg.norm(
                                np.cross(start[a] - start[b], start[c2] - start[b])
                            ),
                            np.linalg.norm(
                                np.cross(start[b] - start[c2], start[d] - start[c2])
                            ),
                        )
                        < 0.3
                    ):
                        continue
                    err = angular_error(
                        angle(start, (a, b, c2, d)), angle(xyz, (a, b, c2, d))
                    )
                    maximum_dihedral = max(maximum_dihedral, err)
                    assert err < 0.5
            for center, neighbors in adjacent.items():
                for subset in itertools.combinations(neighbors, 3):
                    before = np.linalg.det([start[i] - start[center] for i in subset])
                    after = np.linalg.det([xyz[i] - xyz[center] for i in subset])
                    if abs(before) > 0.01:
                        assert before * after > 0
            goal = 180 * round(phi0 / 180) + (180 if form == "far" else 0)
            err = angular_error(angle(xyz, site), goal)
            maximum_target = max(maximum_target, err)
            assert err < 0.12
            for other in adjacent[site[2]] - {site[1], site[3]}:
                other_angle = angle(xyz, (o, c, n, other))
                pyramidal.append(
                    (abs((other_angle + 90) % 180 - 90), path.stem, form, other)
                )
            nstates += 1
    print(
        f"{nstates} locked inputs; {already_fixed} native trees already locked. Max rigid distance error {maximum_geometry:.6f} A, other non-collinear torsion {maximum_dihedral:.5f} deg, selected planar error {maximum_target:.5f} deg; other N substituent greatest departure from plane {max(pyramidal)}",
        flush=True,
    )


def crystal_reference(base, lig, rec, mol, serial, site):
    path = base / "fdv-second"
    source = audited_reference(path, lig)
    perms = source.GetSubstructMatches(mol, uniquify=False, maxMatches=100000)
    assert 0 < len(perms) < 100000
    ids = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[O-]"))
    for _, left, right in ids:
        swapped = []
        for mp in perms:
            ix = list(mp)
            ix[left], ix[right] = ix[right], ix[left]
            swapped.append(tuple(ix))
        perms = (*perms, *swapped)
    src = backbone(path, lig)
    dst = backbone(path, rec)
    keys = sorted(set(src) & set(dst))
    a, b = [np.array([d[i] for i in keys]) for d in (src, dst)]
    rot, _ = Rotation.align_vectors(b - b.mean(axis=0), a - a.mean(axis=0))
    ref = rot.apply(source.GetConformer().GetPositions() - a.mean(axis=0)) + b.mean(
        axis=0
    )
    refs = ref[np.array(perms)]
    ((ni,),) = mol.GetSubstructMatches(Chem.MolFromSmarts("[N+;H1;D3]"))
    ordered = [x.GetIdx() for x in mol.GetAtomWithIdx(ni).GetNeighbors()]

    def sign(arr):
        return np.sign(np.linalg.det(arr[ordered] - arr[ni]))

    assert all(sign(r) == sign(refs[0]) for r in refs)
    return refs, sign, angle(dict(zip(serial, refs[0], strict=True)), site)


def raw_pools(base, matrix, folder):
    data = {}
    for r in csv.DictReader((folder / "modes.csv").open()):
        key = (
            r["panel"],
            r["rec"],
            r["ligand"],
            r["state"],
            r["seed"],
            r["slot"],
            r["form"],
            int(r["mode"]),
        )
        assert key not in data
        data[key] = r
    archived = data.copy()
    refs = {}
    rms_err = phi_err = xtal_err = 0.0
    for path in sorted(matrix.glob("*.pdbqt")):
        panel, rec, lig, state, seed, slot, form = path.stem.split("_")
        log = path.with_suffix(".log").read_text()
        assert (
            f"Performing docking (random seed: {int(seed) + 1009 * int(slot)})" in log
            and log.count("Search trajectory:") == 1
        )
        free = (
            (matrix / "inputs" / f"{panel}_{lig}_{state}_free.pdbqt")
            .read_text()
            .splitlines()
        )
        mol, serial, _, site = geometry(free)
        if panel == "ca":
            key = lig, rec
            if key not in refs:
                refs[key] = crystal_reference(base, lig, rec, mol, serial, site)
        for block in path.read_text().split("ENDMDL")[:-1]:
            lines = block.strip().splitlines()
            mode = int(lines[0].split()[1])
            xyz = atoms(lines)
            r = data.pop((panel, rec, lig, state, seed, slot, form, mode))
            if site:
                phi = angle(xyz, site)
                phi_err = max(phi_err, angular_error(float(r["phi"]), phi))
            else:
                assert not r["phi"]
            if panel == "ca":
                ref, sign, xtal = refs[lig, rec]
                coords = np.array([xyz[i] for i in serial])
                rms = np.sqrt(np.square(ref - coords).sum(axis=2).mean(axis=1)).min()
                rms_err = max(rms_err, abs(float(r["rmsd"]) - rms))
                xtal_err = max(xtal_err, angular_error(float(r["crystal_phi"]), xtal))
                assert int(r["hand"]) == (sign(coords) == sign(ref[0]))
    assert not data and rms_err < 2e-5 and phi_err < 1e-7 and xtal_err < 0.01, (
        rms_err,
        phi_err,
        xtal_err,
    )
    print(
        f"{len(archived)} mode identities/dihedrals, {len(refs)} deposited ligand-receptor alignments; max RMSD / amide / crystal dihedral difference {rms_err:.3g} A / {phi_err:.3g} / {xtal_err:.3g} deg",
        flush=True,
    )
    return archived


def check_cells(modes, folder):
    cells = list(csv.DictReader((folder / "cells.csv").open()))
    groups = defaultdict(list)
    for key, r in modes.items():
        groups[key[:3] + key[4:5]].append(r)
    strata = defaultdict(lambda: [0, 0, 0])
    for cell in cells:
        panel, rec, lig, seed, arm = map(
            cell.get, ("panel", "receptor", "ligand", "seed", "arm")
        )
        group = groups[panel, rec, lig, seed]
        floor = min(
            (r for r in group if r["form"] == "free"),
            key=lambda r: float(r["inter"]) + float(r["intra"]),
        )
        U = float(floor["intra"])
        acceptable = {"both": ("near", "far"), "reverse": ("far", "near")}.get(
            arm, (arm, arm)
        )
        if all(not r["phi"] for r in group):
            acceptable = ("free", "free")
        group = [r for r in group if r["form"] == acceptable[int(r["slot"])]]
        sign = lambda r: int(float(r["rmsd"]) <= 2 and int(r["hand"]) == 1)
        E = lambda r: (
            (float(r["inter"]) + float(r["intra"]) - U) * float(r["slope"])
            + int(cell["barrier"]) * math.sin(math.radians(float(r["phi"] or 0))) ** 2
        )
        ranked = min(group, key=E)
        assert (
            all(ranked[k] == cell[k] for k in ("state", "slot", "mode"))
            and abs(E(ranked) - float(cell["score"])) < 1e-8
        )
        if panel != "fxa":
            assert sign(ranked) == int(cell["strict"])
            poolhit = max(map(sign, group))
            assert poolhit == int(cell["oracle"])
            if panel == "ca" and (arm == "free" or cell["barrier"] == "0"):
                stat = strata[arm, cell["barrier"], rec == lig]
                stat[0] += 1
                stat[1] += sign(ranked)
                stat[2] += poolhit
    for r in csv.DictReader((folder / "summary.csv").open()):
        group = [
            c for c in cells if all(c[k] == r[k] for k in ("panel", "arm", "barrier"))
        ]
        assert len(group) == int(r["cells"])
        for k, l in (("top", "strict"), ("oracle", "oracle")):
            if r[k]:
                assert int(r[k]) == sum(int(c[l]) for c in group)
        gaps = defaultdict(list)
        for a, b in itertools.permutations(group, 2):
            if (
                a["seed"] == b["seed"]
                and a["receptor"] == b["receptor"]
                and float(b["measurement"]) >= 3 * float(a["measurement"])
            ):
                gaps[a["ligand"], b["ligand"]].append(
                    float(b["score"]) - float(a["score"])
                )
        assert len(gaps) == int(r["pairs"])
        assert sum(
            (1 + np.sign(statistics.median(v))) / 2 for v in gaps.values()
        ) == float(r["unique"])
    print(
        f"Verified {len(cells)} top selections, candidate hits and all unique-pair summaries; hCA II (cells,strict,oracle) by (arm,barrier,cognate): {dict(strata)}",
        flush=True,
    )


def verify_energies(binary, base, matrix, records, workers):
    def pool(path):
        panel, rec, lig, state, seed, slot, form = path.stem.split("_")
        if panel == "fxa":
            receptor = (
                base / ("fdv-fxa1" if rec == "1FJS" else "fdv-fxa") / f"{rec}.pdbqt"
            )
        else:
            receptor = (
                base
                / ("fdv-second" if panel == "ca" else "fdv-brd4")
                / "prepared"
                / f"{rec}_receptor.pdbqt"
            )
        inputs = matrix / "inputs"
        ligand = inputs / f"{panel}_{lig}_orig_free.pdbqt"
        bonds = inputs / f"{panel}_{lig}_orig.bonds"
        original = ligand.read_text().splitlines()
        maximum = 0.0
        n = 0
        with tempfile.TemporaryDirectory() as td:
            fixed = Path(td) / "position.pdbqt"
            for block in path.read_text().split("ENDMDL")[:-1]:
                lines = block.strip().splitlines()
                mode = int(lines[0].split()[1])
                r = records[panel, rec, lig, state, seed, slot, form, mode]
                fields = {int(l[6:11]): l[30:54] for l in lines if atom_line(l)}
                fixed.write_text(
                    "\n".join(
                        l[:30] + fields[int(l[6:11])] + l[54:] if atom_line(l) else l
                        for l in original
                    )
                    + "\n"
                )
                out = subprocess.check_output(
                    [
                        binary,
                        "--cpu",
                        "1",
                        "--verbosity",
                        "0",
                        "--score_only",
                        "--autobox",
                        "--unbound_energy",
                        "0",
                        "--ligand",
                        ligand,
                        "--ligand_bonds",
                        bonds,
                        "--receptor",
                        receptor,
                        "--pose_coordinates",
                        fixed,
                    ],
                    text=True,
                )
                read = lambda name: float(
                    re.search(r"^\s*" + name + r"\s*:\s*([-\d.]+)", out, re.M)[1]
                )
                inter, intra, slope = [float(r[c]) for c in ("inter", "intra", "slope")]
                err = max(
                    abs(read("Ligand - Receptor") - inter),
                    abs(read("Ligand") - intra),
                    abs(
                        read("Estimated Free Energy of Binding")
                        - (inter + intra) * slope
                    ),
                )
                assert err <= 0.000501, (path, mode, err)
                maximum = max(maximum, err)
                n += 1
        return n, maximum

    with ThreadPoolExecutor(workers) as executor:
        totals = list(executor.map(pool, sorted(matrix.glob("*.pdbqt"))))
    print(
        f"Standard autoboxed --score_only, explicit zero: {sum(n for n, _ in totals)} independently rescored rounded structures, worst INTER/INTRA/standard-score discrepancy {max(e for _, e in totals):.8f} kcal/mol (ordinary CLI prints three decimals).",
        flush=True,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("base", "matrix", "results"):
        p.add_argument(arg, type=Path)
    p.add_argument("--inputs-only", action="store_true")
    p.add_argument(
        "--binary",
        type=Path,
        help="also recompute all common-topology scores through ordinary score_only",
    )
    p.add_argument("--jobs", type=int, default=3)
    a = p.parse_args()
    audit_inputs(a.base, a.matrix)
    if not a.inputs_only:
        rows = raw_pools(a.base, a.matrix, a.results)
        check_cells(rows, a.results)
        if a.binary:
            verify_energies(a.binary.resolve(), a.base, a.matrix, rows, a.jobs)


if __name__ == "__main__":
    main()
