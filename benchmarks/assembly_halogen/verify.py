"""Independent replay through Meeko, receptor PDBQT records, and scalar geometry."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np
from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate


def pxyz(line):
    return tuple(map(float, (line[30:38], line[38:46], line[46:54])))


def oxygen_atoms(folder):
    receptor = (folder / "receptor.pdbqt").read_text()
    accepted = []
    terminal = [
        q
        for q in receptor.splitlines()
        if q.startswith("ATOM") and q[12:16].strip() == "OXT"
    ]
    for l in receptor.splitlines():
        if l.startswith("ATOM") and l[77:].strip() == "OA":
            residue = l[17:20].strip()
            atom = l[12:16].strip()
            if (
                (
                    atom == "O"
                    and not any(
                        q[17:27] == l[17:27] and q[12:16].strip() == "OXT"
                        for q in terminal
                    )
                )
                or (residue == "ASN" and atom == "OD1")
                or (residue == "GLN" and atom == "OE1")
            ):
                accepted.append(pxyz(l))
    return accepted


def scalar_scores(mol, conf, acceptors):
    accum = [0.0, 0.0]
    for atom in mol.GetAtoms():
        elem = atom.GetAtomicNum()
        if elem not in (17, 35, 53):
            continue
        n = atom.GetNeighbors()
        if len(n) != 1 or n[0].GetAtomicNum() != 6:
            continue
        r0, eps = {17: (3.27, 0.265), 35: (3.37, 0.32), 53: (3.50, 0.40)}[elem]
        x = tuple(conf.GetAtomPosition(atom.GetIdx()))
        c = tuple(conf.GetAtomPosition(n[0].GetIdx()))
        best = [0.0, 0.0]
        for o in acceptors:
            u = np.subtract(c, x)
            v = np.subtract(o, x)
            d = math.dist(x, o)
            if abs(d - r0) > 4:
                continue
            theta = math.degrees(
                math.atan2(np.linalg.norm(np.cross(u, v)), np.dot(u, v))
            )
            z = min(1, max(0, (theta - 135) / 30))
            g = math.exp(-((d - r0) ** 2) / (2 * 0.35**2))
            best[0] = max(best[0], g)
            best[1] = max(best[1], g * (3 * z * z - 2 * z**3))
        accum = [t - eps * b for t, b in zip(accum, best)]
    return accum


def graph_maps(probe, ref, chirality):
    matches = ref.GetSubstructMatches(
        probe, uniquify=False, useChirality=chirality, maxMatches=100000
    )
    if not matches:
        raise ValueError("no native atom correspondence")
    # CalcRMS exchanges terminal O atoms in carboxyl/nitro groups.
    swaps = []
    for center in ref.GetAtoms():
        if center.GetAtomicNum() not in (6, 7):
            continue
        oxy = [
            n.GetIdx()
            for n in center.GetNeighbors()
            if n.GetAtomicNum() == 8 and n.GetDegree() == 1
        ]
        if len(oxy) == 2:
            swaps.append(oxy)
    all_maps = set(matches)
    for a, b in swaps:
        for mapping in list(all_maps):
            all_maps.add(tuple(b if i == a else a if i == b else i for i in mapping))
    return list(all_maps)


def replay(folder, pool, binary, out):
    check = out / (folder.name + "-" + pool.name)
    seed = pool.stem.split("-")[-1]
    subprocess.run(
        [
            binary,
            "--receptor",
            folder / "receptor.pdbqt",
            "--ligand",
            folder / "query.pdbqt",
            "--config",
            folder / "box",
            "--seed",
            seed,
            "--cpu",
            "1",
            "--exhaustiveness",
            "3",
            "--max_evals",
            "60000",
            "--num_modes",
            "9",
            "--pose_selection",
            "late",
            "--refine_candidates",
            "12",
            "--energy_range",
            "100",
            "--out",
            check,
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if pool.read_bytes() != check.read_bytes():
        raise ValueError(f"changed default replay for {pool}")
    return hashlib.sha256(check.read_bytes()).hexdigest()


def main(work, results, binary, replayed):
    with (results / "poses.csv").open() as f:
        expected = {
            (r["case"], int(r["seed"]), int(r["model"])): r for r in csv.DictReader(f)
        }
    n = 0
    rms_error = 0.0
    energy_error = 0.0
    pools = []
    chiral = []
    for folder in sorted(work.iterdir()):
        if not (folder / "ready").exists() or not list(folder.glob("poses-*.pdbqt")):
            continue
        ref = Chem.SDMolSupplier(str(folder / "crystal.sdf"), removeHs=True)[0]
        acceptors = oxygen_atoms(folder)
        for pool in sorted(folder.glob("poses-*.pdbqt")):
            pools.append((folder, pool))
            raw = PDBQTMolecule.from_file(str(pool), skip_typing=True)
            mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(raw)[0])
            maps = np.array(graph_maps(mol, ref, False))
            stereo_refs = ref.GetConformer().GetPositions()[
                np.array(graph_maps(mol, ref, True))
            ]
            refs = ref.GetConformer().GetPositions()[maps]
            scores = [
                float(l.split()[3])
                for l in pool.read_text().splitlines()
                if l.startswith("REMARK VINA RESULT:")
            ]
            for i, c in enumerate(mol.GetConformers(), 1):
                row = expected[(folder.name, int(pool.stem.split("-")[-1]), i)]
                rms = np.sqrt(
                    np.square(refs - c.GetPositions()).sum(axis=2).mean(axis=1)
                ).min()
                rms_error = max(rms_error, abs(float(row["rmsd"]) - rms))
                stereo = np.sqrt(
                    np.square(stereo_refs - c.GetPositions()).sum(axis=2).mean(axis=1)
                ).min()
                chiral.append(
                    dict(
                        case=folder.name,
                        seed=int(pool.stem.split("-")[-1]),
                        model=i,
                        stereo_rmsd_A=stereo,
                        stereo_delta_A=stereo - rms,
                        threshold_changed=int((stereo <= 2) != (rms <= 2)),
                    )
                )
                if scores[i - 1] != float(row["native"]):
                    raise ValueError("native score differs")
                for v, k in zip(
                    scalar_scores(mol, c, acceptors),
                    ("radial_delta", "directional_delta"),
                ):
                    energy_error = max(energy_error, abs(v - float(row[k])))
                n += 1
    if n != len(expected) or rms_error > 2e-5 or energy_error > 2e-5:
        raise ValueError(
            f"{n} models; max RMSD disagreement {rms_error}, energy disagreement {energy_error}"
        )
    with (results / "stereo_sensitivity.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=chiral[0])
        w.writeheader()
        w.writerows(chiral)
    output = dict(
        models=n,
        pools=len(pools),
        max_RMSD_difference_A=rms_error,
        max_correction_difference=energy_error,
        max_stereo_map_difference_A=max(v["stereo_delta_A"] for v in chiral),
        models_changing_stereo_threshold=sum(v["threshold_changed"] for v in chiral),
    )
    if binary:
        replayed.mkdir(exist_ok=True, parents=True)
        with ThreadPoolExecutor(3) as ex:
            hashes = list(ex.map(lambda pair: replay(*pair, binary, replayed), pools))
        output["byte_identical_native_replays"] = len(hashes)
        output["replay_binary_sha256"] = hashlib.sha256(binary.read_bytes()).hexdigest()
    (results / "audit.json").write_text(json.dumps(output, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("results", type=Path)
    p.add_argument("--binary", type=Path)
    p.add_argument("--replayed", type=Path)
    a = p.parse_args()
    main(a.work, a.results, a.binary.resolve() if a.binary else None, a.replayed)
