"""Independent PDB SMTRY/CRYST1 expansion and RCSB-instance RDKit checks of packing and RMSD."""

import argparse
import csv
import itertools
import math
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from scipy.spatial import cKDTree


def lines(path):
    return path.read_text().splitlines()


def pos(line):
    return [float(line[k : k + 8]) for k in (30, 38, 46)]


def images(path, host, probes):
    raw = lines(path)
    cell = next(l for l in raw if l.startswith("CRYST1"))
    a, b, c = [float(cell[i:j]) for i, j in ((6, 15), (15, 24), (24, 33))]
    alpha, beta, gamma = np.deg2rad(
        [float(cell[i:j]) for i, j in ((33, 40), (40, 47), (47, 54))]
    )
    cy = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    axes = np.array(
        [
            [a, 0, 0],
            [b * np.cos(gamma), b * np.sin(gamma), 0],
            [c * np.cos(beta), c * cy, c * np.sqrt(1 - np.cos(beta) ** 2 - cy**2)],
        ]
    )
    ops = defaultdict(list)
    for line in raw:
        if line.startswith("REMARK 290   SMTRY"):
            _, _, _, num, *x = line.split()
            ops[int(num)].append(list(map(float, x)))
    rows = [
        r
        for r in raw
        if r.startswith("ATOM  ")
        and r[76:78].strip() not in ("H", "D")
        and float(r[54:60]) > 0
    ]
    xyz = np.array([pos(r) for r in rows])
    other = np.array([r[21] != host for r in rows])
    center = np.mean(probes, axis=0)
    radius = np.max(np.linalg.norm(probes - center, axis=1)) + 20
    xyz_parts = []
    atoms = []
    for op, matrix in ops.items():
        mat = np.asarray(matrix)
        trans = xyz @ mat[:, :3].T + mat[:, 3]
        for n in itertools.product(range(-3, 4), repeat=3):
            shift = np.array(n) @ axes
            if (
                np.linalg.norm(center - trans.mean(0) - shift)
                > radius + np.linalg.norm(trans - trans.mean(0), axis=1).max()
            ):
                continue
            x = trans + shift
            keep = np.linalg.norm(x - center, axis=1) < radius
            if op == 1 and n == (0, 0, 0):
                keep &= other
            xyz_parts.append(x[keep])
            atoms.extend(
                (rows[i][76:78].strip(), float(rows[i][54:60]))
                for i in np.flatnonzero(keep)
            )
    return np.vstack(xyz_parts), atoms


def coordinate_sets(path):
    sets = []
    coords = []
    elements = []
    for l in lines(path):
        if l.startswith("MODEL"):
            coords = []
            elements = []
        if l.startswith(("ATOM", "HETATM")):
            kind = l.split()[-1]
            if kind not in ("H", "HD", "HS"):
                coords.append(pos(l))
                elements.append(
                    {"A": "C", "OA": "O", "NA": "N", "SA": "S"}.get(kind, kind.upper())
                )
        if l.startswith("ENDMDL"):
            sets.append((np.array(coords), elements))
    return sets


def protein_alignment(source, schain, target, tchain):
    def ca(path, chain):
        return {
            int(l[22:26]): pos(l)
            for l in lines(path)
            if l.startswith("ATOM")
            and l[21] == chain
            and l[12:16].strip() == "CA"
            and l[16] in " A"
            and float(l[54:60]) > 0
            and 44 <= int(l[22:26]) <= 166
        }

    sa, tb = ca(source, schain), ca(target, tchain)
    ids = sorted(sa.keys() & tb.keys())
    a = np.array([sa[k] for k in ids])
    b = np.array([tb[k] for k in ids])
    u, _, vt = np.linalg.svd((a - a.mean(0)).T @ (b - b.mean(0)))
    sign = np.linalg.det(u @ vt)
    return u @ np.diag([1, 1, sign]) @ vt, a.mean(0), b.mean(0)


def molecule(path, charged):
    mol = Chem.SDMolSupplier(str(path), removeHs=True)[0]
    if charged:
        (q,) = mol.GetSubstructMatches(Chem.MolFromSmarts("[N;X3;R;!$(N-C=O);!$(N-c)]"))
        at = mol.GetAtomWithIdx(q[0])
        at.SetFormalCharge(1)
        at.SetNumExplicitHs(1)
        Chem.SanitizeMol(mol)
    return mol


def main(a):
    with (a.out / "candidates.csv").open() as f:
        records = list(csv.DictReader(f))
    p_by_r = defaultdict(list)
    for r in records:
        p_by_r[r["target"], r["receptor"]].append(r)
    worst_r = worst_d = worst_e = 0
    n_r = n_e = 0
    vals = dict(C=1.9, N=1.8, O=1.7, S=2.0, CL=1.8, F=1.5)
    for (target, receptor), rs in p_by_r.items():
        files = defaultdict(list)
        for r in rs:
            files[r["ligand"], r["seed"], r["slot"]].append(r)
        path = (a.ca if target == "hcaii" else a.brd) / "raw" / (receptor + ".pdb")
        inputdir = a.ca_runs if target == "hcaii" else a.brd_runs
        chain = "B" if receptor == "4OGJ" else "A"
        pool = {
            (l, s, t): coordinate_sets(inputdir / f"{receptor}_{l}_{s}_{t}.pdbqt")
            for l, s, t in files
        }
        candidates = np.vstack([x for mols in pool.values() for x, _ in mols])
        xyz, mates = images(path, chain, candidates)
        tree = cKDTree(xyz)
        for (lig, s, t), subset in files.items():
            if target == "brd4":
                source = a.brd / "raw" / (lig + ".pdb")
                copies = []
                comp = {"4OGI": "R78", "4OGJ": "2TA", "4O76": "1M3", "4O77": "2RE"}[lig]
                ids = sorted(
                    {
                        (l[21], l[22:26].strip())
                        for l in lines(source)
                        if l.startswith("HETATM") and l[17:20] == comp
                    }
                )
                for schain, num in ids:
                    local = a.instances / f"{lig}_{schain}_{num}.sdf"
                    if not local.exists():
                        url = f"https://models.rcsb.org/v1/{lig}/ligand?auth_asym_id={schain}&auth_seq_id={num}&encoding=sdf"
                        local.write_bytes(urlopen(url, timeout=60).read())
                    mol = molecule(local, lig != "4O77")
                    rot, src, dst = protein_alignment(source, schain, path, chain)
                    conf = mol.GetConformer()
                    conf.SetPositions((conf.GetPositions() - src) @ rot + dst)
                    copies.append(mol)
                file = inputdir / f"{receptor}_{lig}_{s}_{t}.pdbqt"
                q = Chem.RemoveHs(
                    RDKitMolCreate.from_pdbqt_mol(
                        PDBQTMolecule.from_file(str(file), skip_typing=True)
                    )[0]
                )
            for r in subset:
                model = int(r["model"]) - 1
                x, elems = pool[lig, s, t][model]
                d, i = tree.query(x)
                worst_d = max(worst_d, abs(float(r["closest"]) - min(d)))
                # Sum each ligand atom's actual neighbors. The 4.5 A upper bound exceeds every radius sum used here.
                energy = 0.0
                for coord, e, neighbors in zip(
                    x, elems, tree.query_ball_point(x, 4.5), strict=True
                ):
                    for j in neighbors:
                        elem, occ = mates[j]
                        energy += (
                            0.840245
                            * occ
                            * max(0.0, vals[e] + vals[elem] - math.dist(coord, xyz[j]))
                            ** 2
                        )
                worst_e = max(worst_e, abs(energy - float(r["repulsion"])))
                n_e += 1
                if target == "brd4":
                    calc = min(
                        rdMolAlign.CalcRMS(q, mol, prbId=model) for mol in copies
                    )
                    worst_r = max(worst_r, abs(calc - float(r["alt_rmsd"])))
                    n_r += 1
    if max(worst_r, worst_d, worst_e) > 2e-4:
        raise ValueError((worst_r, worst_d, worst_e))
    print(
        f"SMTRY and explicit unit-cell translations reproduced {n_e} candidate nearest distances and repulsion energies to {worst_d:.3g} A and {worst_e:.3g} raw units"
    )
    print(
        f"RDKit and independently obtained ligand-instance SDFs reproduced {n_r} best-copy RMSDs to {worst_r:.3g} A"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("ca", "brd", "ca_runs", "brd_runs", "out", "instances"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    a.instances.mkdir(exist_ok=True)
    main(a)
