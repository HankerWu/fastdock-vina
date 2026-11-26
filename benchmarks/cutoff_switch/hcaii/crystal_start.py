"""Project cognate X-ray ligands onto the prepared torsion tree and locally refine them."""

import argparse
import csv
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from prep import CASES, RECS, xyz
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
from score_physical_pool import mapped_ligand_bonds


def atom_coordinates(text):
    return {
        int(line[6:11]): xyz(line)
        for line in text.splitlines()
        if line.startswith(("ATOM", "HETATM"))
    }


def torsion_tree(text):
    torsions, open_branches = [], []
    for line in text.splitlines():
        if line.startswith("BRANCH "):
            _, first, second = line.split()
            branch = (int(first), int(second), [])
            torsions.append(branch)
            open_branches.append(branch)
        elif line.startswith("ENDBRANCH "):
            open_branches.pop()
        elif line.startswith(("ATOM", "HETATM")):
            for branch in open_branches:
                branch[2].append(int(line[6:11]))
    return torsions


def align_to_reference(coords, heavy, target):
    center = coords[heavy].mean(axis=0)
    target_center = target.mean(axis=0)
    left, _, right = np.linalg.svd(
        (coords[heavy] - center).T @ (target - target_center)
    )
    left[:, -1] *= np.sign(np.linalg.det(left @ right))
    return (coords - center) @ (left @ right) + target_center


def project(prepared, serials, reference, edges):
    """Only rigid placement and BRANCH rotations are fitted. Fragment geometry stays fixed."""
    coords = atom_coordinates(prepared)
    initial = np.zeros((max(coords) + 1, 3))
    for serial, point in coords.items():
        initial[serial] = point
    tree = torsion_tree(prepared)
    heavy = set(serials)
    target = dict(zip(serials, reference))
    adjacency = {serial: set() for serial in coords}
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)

    def frame_angle(a, u, v, d):
        axis = v - u
        axis /= np.linalg.norm(axis)
        x = a - u
        y = d - v
        x -= np.dot(axis, x) * axis
        y -= np.dot(axis, y) * axis
        return np.arctan2(np.dot(axis, np.cross(x, y)), np.dot(x, y))

    guess = np.zeros(len(tree))
    for i, (u, v, mobile) in enumerate(tree):
        left = sorted((adjacency[u] & heavy) - {v})
        right = sorted((adjacency[v] & heavy) - {u})
        if left and right:
            order = [left[0], u, v, right[0]]
            guess[i] = frame_angle(*(target[n].copy() for n in order)) - frame_angle(
                *(coords[n].copy() for n in order)
            )

    def conformation(angles):
        moved = initial.copy()
        for theta, (u, v, indices) in zip(angles, tree):
            axis = moved[v] - moved[u]
            rotation = Rotation.from_rotvec(theta * axis / np.linalg.norm(axis))
            moved[indices] = rotation.apply(moved[indices] - moved[u]) + moved[u]
        return align_to_reference(moved, serials, reference)

    fit = least_squares(
        lambda angles: (conformation(angles)[serials] - reference).ravel(),
        guess,
        max_nfev=250,
    )
    if not fit.success:
        raise ValueError(f"Torsion projection failed: {fit.message}")
    placed = conformation(fit.x)
    lines = []
    for line in prepared.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            pt = placed[int(line[6:11])]
            line = line[:30] + "".join(f"{v:8.3f}" for v in pt) + line[54:]
        lines.append(line)
    return "\n".join(lines) + "\n"


def energy(common, base, pose, flag):
    text = subprocess.check_output(
        common
        + [
            "--ligand",
            str(base),
            "--pose_coordinates",
            str(pose),
            "--score_only",
            "--net_polar",
            "native",
            "--verbosity",
            "0",
        ]
        + flag,
        text=True,
    )
    inter, intra, _ = map(
        float, re.search(r"^Net polar native: (.+)$", text, re.M)[1].split()
    )
    return inter, intra


def run_ligand(binary, work, output, ligand):
    prepared = work / "prepared"
    meta = json.loads((prepared / "metadata.json").read_text())[ligand]
    reference = np.array(meta["references"][ligand]["coords"])
    serials = meta["qt_serial"]
    base = prepared / f"{ligand}_ligand.pdbqt"
    edges = mapped_ligand_bonds(base.read_text())
    bonds = output / f"{ligand}.bonds"
    bonds.write_text("".join(f"{a} {b}\n" for a, b in edges))
    seedpose = output / f"{ligand}_start.pdbqt"
    seedpose.write_text(project(base.read_text(), serials, reference, edges))

    ref = Chem.SDMolSupplier(str(prepared / f"{ligand}_graph.sdf"))[0]
    for i, pt in enumerate(reference):
        ref.GetConformer().SetAtomPosition(i, pt)

    def rms(path):
        coords = atom_coordinates(path.read_text())
        pr = Chem.Mol(ref)
        for i, serial in enumerate(serials):
            pr.GetConformer().SetAtomPosition(i, coords[serial])
        return rdMolAlign.CalcRMS(
            pr, ref, maxMatches=1000000, symmetrizeConjugatedTerminalGroups=True
        )

    common = [
        str(binary),
        "--receptor",
        str(prepared / f"{ligand}_receptor.pdbqt"),
        "--config",
        str(prepared / f"{ligand}_box.txt"),
        "--ligand_bonds",
        str(bonds),
        "--cpu",
        "1",
        "--seed",
        "42",
    ]
    rows = []
    for arm, flag in [("native", []), ("smooth", ["--smooth_cutoff"])]:
        refined = output / f"{ligand}_{arm}_refined.pdbqt"
        subprocess.run(
            common
            + [
                "--ligand",
                str(seedpose),
                "--local_only",
                "--out",
                str(refined),
                "--verbosity",
                "0",
            ]
            + flag,
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        )
        projected_energy = energy(common, base, seedpose, flag)
        refined_energy = energy(common, base, refined, flag)
        for seed in [42, 1729]:
            pool_dir = work / ("runs" if ligand in RECS else "cognate")
            pool = (pool_dir / f"{ligand}_{ligand}_{seed}_native.pdbqt").read_text()
            modes = re.findall(r"^MODEL\s+\d+\s*\n.*?^ENDMDL$", pool, re.M | re.S)
            assert len(modes) > 0
            pool_values = []
            for m, block in enumerate(modes, 1):
                pose = output / f"{ligand}_{seed}_{m}.pdbqt"
                pose.write_text(block + "\n")
                pool_values.append(energy(common, base, pose, flag))
            best_index = min(range(len(pool_values)), key=lambda i: sum(pool_values[i]))
            pool_energy = pool_values[best_index]
            rows.append(
                dict(
                    ligand=ligand,
                    arm=arm,
                    seed=seed,
                    projected_RMSD=round(rms(seedpose), 6),
                    refined_RMSD=round(rms(refined), 6),
                    projected_inter_intra=round(sum(projected_energy), 6),
                    refined_inter_intra=round(sum(refined_energy), 6),
                    n_modes=len(modes),
                    best_pool_model=best_index + 1,
                    best_pool_inter_intra=round(sum(pool_energy), 6),
                    refined_minus_best=round(sum(refined_energy) - sum(pool_energy), 6),
                    delta_inter=round(refined_energy[0] - pool_energy[0], 6),
                    delta_intra=round(refined_energy[1] - pool_energy[1], 6),
                )
            )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    output = args.work / "crystal_start"
    output.mkdir()
    with ThreadPoolExecutor(args.jobs) as pool:
        futures = [
            pool.submit(run_ligand, args.binary.resolve(), args.work, output, lig)
            for lig in CASES
        ]
        rows = [r for f in futures for r in f.result()]
    with (output / "results.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
