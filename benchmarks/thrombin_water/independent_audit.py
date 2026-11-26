"""Recompute RMSD and both scores directly from raw PDBQT using wwPDB CCD bond tables."""

import csv
import itertools
import json
import math
from collections import Counter
from pathlib import Path
from statistics import median

import gemmi
import numpy as np
from rdkit import Chem
from scipy.spatial.transform import Rotation


# Deliberately no imports from the experiment/evaluation scripts.
def read_csv(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def positions(line):
    return [float(line[30:38]), float(line[38:46]), float(line[46:54])]


def raw_lines(pdb):
    return Path("raw", pdb + ".pdb").read_text().splitlines()


def landmarks(pdb):
    ids = ("57", "60A", "60D", "99", "174", "189", "195", "215", "216", "227", "228")
    selected = {
        l[22:27].strip(): positions(l)
        for l in raw_lines(pdb)
        if l.startswith("ATOM")
        and l[12:16].strip() == "CA"
        and l[21] == "H"
        and l[16] in " A"
    }
    return np.array([selected[i] for i in ids])


def move_points(coords, source, receiver):
    a, b = landmarks(source), landmarks(receiver)
    rotation, _ = Rotation.align_vectors(b - b.mean(0), a - a.mean(0))
    return rotation.apply(coords - a.mean(0)) + b.mean(0)


def refs_from_ccd(lig_pdb, component, topology, receiver):
    # CCD named bonds provide adjacency without inferring bonds from measured coordinates.
    pdb = {
        r[12:16].strip(): positions(r)
        for r in raw_lines(lig_pdb)
        if r.startswith("HETATM") and r[17:20] == component and r[16] in " A"
    }
    cif = gemmi.cif.read_file("raw/" + component + ".cif").sole_block()
    ccd = Chem.RWMol()
    idx = {}
    for name, element in cif.find("_chem_comp_atom.", ["atom_id", "type_symbol"]):
        if name not in pdb:
            continue
        idx[name] = ccd.AddAtom(Chem.Atom(element.title()))
    for a, b in cif.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2"]):
        if a in idx and b in idx:
            ccd.AddBond(idx[a], idx[b], Chem.BondType.SINGLE)
    # Canonical SMILES and SMILES IDX use a different atom order from the original input.
    canonical = Chem.MolFromSmiles(
        next(l[14:] for l in topology if l.startswith("REMARK SMILES "))
    )
    flat = Chem.RWMol()
    for a in canonical.GetAtoms():
        flat.AddAtom(Chem.Atom(a.GetAtomicNum()))
    for b in canonical.GetBonds():
        flat.AddBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), Chem.BondType.SINGLE)
    match = ccd.GetSubstructMatches(flat, uniquify=False)
    if (
        not match
        or ccd.GetNumAtoms() != flat.GetNumAtoms()
        or ccd.GetNumBonds() != flat.GetNumBonds()
    ):
        raise ValueError(component)
    raw = np.array([pdb[name] for name in idx])
    return np.array([move_points(raw[list(m)], lig_pdb, receiver) for m in match])


def get_serials(prologue):
    idx = {}
    for r in prologue:
        if r.startswith("REMARK SMILES IDX"):
            val = [int(x) for x in r.split()[3:]]
            idx.update(zip(val[::2], val[1::2]))
    return [idx[i] for i in range(1, len(idx) + 1)]


def main():
    detail = {
        (r["rec"], r["ligand"], r["seed"], r["model"]): r for r in read_csv("poses.csv")
    }
    series = {r["id"]: r for r in read_csv("series.csv")}
    ref_cache = {}
    topology = {}
    pool_best = {}
    pose_success = Counter()
    pose_non = Counter()
    max_rms = max_energy = max_water = 0.0
    models = struct = 0
    site = [
        positions(l)
        for l in raw_lines("3D49")
        if l.startswith("HETATM") and l[17:26] == "HOH H2006"
    ]
    if len(site) != 1:
        raise ValueError("Expected H:2006")
    for path in sorted(Path("results").glob("*.pdbqt")):
        rec, lig, seed = path.stem.split("_")
        if lig not in topology:
            topology[lig] = Path("prepared", lig + ".pdbqt").read_text().splitlines()
        serials = get_serials(topology[lig])
        p = series[lig]["PDB"]
        cc = series[lig]["CCD"]
        if p and (lig, rec) not in ref_cache:
            ref_cache[lig, rec] = refs_from_ccd(p, cc, topology[lig], rec)
        water = move_points(site, "3D49", rec)[0]
        candidate = []
        for block in path.read_text().split("ENDMDL"):
            vals = block.splitlines()
            try:
                idx = int(next(l for l in vals if l.startswith("MODEL")).split()[1])
            except StopIteration:
                continue
            e = float(
                next(l for l in vals if l.startswith("REMARK VINA RESULT")).split()[3]
            )
            serial_xyz = {
                int(l[6:11]): positions(l)
                for l in vals
                if l.startswith(("ATOM", "HETATM"))
            }
            pos = np.array([serial_xyz[k] for k in serials])
            observed = min(math.dist(v, water) for v in pos)
            t = min(1, max(0, (observed - 2.0) / 0.6))
            corrected = e - (0.00831446261815324 * 298.15 / 4.184) * math.log(
                1 - 0.75 * (1 - (3 * t * t - 2 * t**3))
            )
            saved = detail[rec, lig, seed, str(idx)]
            max_water = max(max_water, abs(float(saved["distance"]) - observed))
            max_energy = max(max_energy, abs(float(saved["exclusion"]) - corrected))
            rms = None
            if p:
                rms = min(
                    math.sqrt(
                        sum(math.dist(u, v) ** 2 for u, v in zip(pos, ref)) / len(pos)
                    )
                    for ref in ref_cache[lig, rec]
                )
                max_rms = max(max_rms, abs(float(saved["rmsd"]) - rms))
                struct += 1
            candidate.append((e, corrected, rms, idx))
            models += 1
        for arm, k in [("native", 0), ("exclusion", 1)]:
            best = min(candidate, key=lambda a: (a[k], a[3]))
            pool_best[arm, rec, seed, lig] = best[k]
            if p:
                pose_success[arm] += best[2] <= 2
                if p != rec:
                    pose_non[arm] += best[2] <= 2
    pairs = {}
    recs = ("3D49", "2ZFF", "2ZC9")
    seeds = ("83519", "194027")
    for arm in ("native", "exclusion"):
        scores = []
        for a, b in itertools.combinations(series, 2):
            a, b = sorted(
                (a, b), key=lambda i: float(series[i]["Kd_uM_derived_at_298.15K"])
            )
            if float(series[b]["Kd_uM_derived_at_298.15K"]) < 3 * float(
                series[a]["Kd_uM_derived_at_298.15K"]
            ):
                continue
            pred = median(
                pool_best[arm, r, s, b] - pool_best[arm, r, s, a]
                for r, s in itertools.product(recs, seeds)
            )
            scores.append(1 if pred > 1e-9 else 0 if pred < -1e-9 else 0.5)
        pairs[arm] = [sum(scores), len(scores)]
    assert len(pool_best) == 144 and len(detail) == models and struct > 0
    assert max(max_rms, max_energy, max_water) < 1e-8
    out = {
        "pools": len(pool_best) // 2,
        "candidate_poses": models,
        "cocrystal_candidate_poses": struct,
        "largest_RMSD_difference_A": max_rms,
        "largest_distance_difference_A": max_water,
        "largest_score_difference_kcal": max_energy,
        "top1_all": dict(pose_success),
        "top1_noncognate": dict(pose_non),
        "unique_Kd_pairs": pairs,
    }
    print(json.dumps(out, indent=2))
    Path("audit.json").write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
