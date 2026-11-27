"""Verify deposited glucose stereochemistry, pose RMSDs and acceptor geometry."""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from acceptor import lp_gate


def explicit_tetrahedron():
    # Four vertices of a regular tetrahedron. Treat the first as O-C and
    # enumerate all permutations of the other three arms as H / LP / LP.
    from itertools import permutations

    vec = np.array(
        [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float
    ) / np.sqrt(3)
    rng = np.random.default_rng(49631)
    rays = rng.standard_normal((64, 3))
    rays /= np.linalg.norm(rays, axis=1)[:, None]
    worst = 0.0
    for h, lp1, lp2 in permutations(vec[1:]):
        for ray in np.concatenate((rays, vec, -vec)):
            actual = lp_gate(np.zeros(3), vec[0], h[None, :], ray)[0]
            exact = max(0.0, ray @ lp1, ray @ lp2) ** 2
            worst = max(worst, abs(actual - exact))
    if worst > 1.0e-12:
        raise ValueError(worst)
    return worst


def xyz(line):
    return np.array([float(line[s : s + 8]) for s in (30, 38, 46)])


def independent_mapping(work, entry, comp):
    raw = (work / "raw" / f"{entry}.pdb").read_text().splitlines()
    query = (work / entry / f"{comp}.pdbqt").read_text().splitlines()
    mol = Chem.MolFromSmiles(
        next(
            l.split()[2]
            for l in query
            if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX")
        )
    )
    source = [
        l[:16] + " " + l[17:]
        for l in raw
        if l.startswith("HETATM")
        and l[21] == "A"
        and l[17:20] == comp
        and l[76:78].strip() not in ["H", "D"]
    ]
    bound = Chem.MolFromPDBBlock(
        "\n".join(source) + "\nEND\n", removeHs=False, sanitize=True
    )
    Chem.AssignStereochemistryFrom3D(bound)
    if Chem.MolToSmiles(bound) != Chem.MolToSmiles(mol):
        raise ValueError(
            f"PDB geometry disagrees with canonical graph for {entry} {comp}"
        )
    serials = {}
    for l in query:
        if l.startswith("REMARK SMILES IDX "):
            seq = list(map(int, l.split()[3:]))
            serials.update(zip(seq[::2], seq[1::2]))
    ids = [serials[i + 1] for i in range(mol.GetNumAtoms())]
    matches = bound.GetSubstructMatches(mol, useChirality=True, uniquify=False)
    refs = np.array([bound.GetConformer().GetPositions()[list(m)] for m in matches])
    if len(refs) != 1:
        raise ValueError(f"Expected a uniquely mapped glucose graph, got {len(refs)}")
    return ids, refs[0]


def check(work, poses):
    with poses.open() as f:
        reported = {
            (r["case"], r["cap"], r["seed"], r["engine"], int(r["rank"])): float(
                r["rmsd"]
            )
            for r in csv.DictReader(f)
        }
    mp = {}
    worst = 0.0
    count = 0
    for pool in sorted((work / "runs").glob("*.pdbqt")):
        case, cap, seed, arm = pool.stem.split("_")
        if case not in mp:
            mp[case] = independent_mapping(work, *case.split("-"))
        ids, ref = mp[case]
        for block in pool.read_text().split("ENDMDL"):
            if not block.strip():
                continue
            ls = block.splitlines()
            rank = int(next(l for l in ls if l.startswith("MODEL")).split()[1])
            px = {int(l[6:11]): xyz(l) for l in ls if l.startswith(("ATOM", "HETATM"))}
            dist = np.linalg.norm(ref - np.array([px[s] for s in ids])) / np.sqrt(
                len(ids)
            )
            worst = max(worst, abs(reported.pop((case, cap, seed, arm, rank)) - dist))
            count += 1
    if reported or worst > 1e-8:
        raise ValueError((len(reported), worst))
    print(
        json.dumps(
            dict(
                verified_poses=count,
                stereochemical_conditions=len(mp),
                max_RMSD_error_A=worst,
                max_tetrahedron_gate_error=explicit_tetrahedron(),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument("--poses", type=Path, default=HERE / "results/poses.csv")
    a = ap.parse_args()
    check(a.work, a.poses)
