"""Measure poses in the receptor frame; check every RMSD against the bound RCSB SDF."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

from external_run import write


def poses(path):
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            result[int(line[6:11])] = (
                line.split()[-1],
                np.asarray([float(line[i : i + 8]) for i in (30, 38, 46)]),
            )
        if line.startswith("ENDMDL"):
            yield result
            result = {}
    if result:
        yield result


def reference(work, pdb):
    d = work / pdb
    meta = json.loads((d / "preparation.json").read_text())
    lines = (d / "ligand.pdbqt").read_text().splitlines()
    order = {}
    for line in lines:
        if line.startswith("REMARK SMILES IDX "):
            ix = list(map(int, line.split()[3:]))
            order.update(zip(ix[::2], ix[1::2]))
    smiles = next(
        l.split()[2]
        for l in lines
        if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX ")
    )
    query = Chem.MolFromSmiles(smiles)
    bound = Chem.MolFromMolFile(str(work / "raw" / (pdb + "-bound.sdf")), removeHs=True)
    if bound.GetNumAtoms() != query.GetNumAtoms() or sorted(order) != list(
        range(1, query.GetNumAtoms() + 1)
    ):
        raise ValueError("Incomplete independent SDF/SMILES atom map")
    Chem.AssignStereochemistryFrom3D(bound)
    matches = bound.GetSubstructMatches(
        query, uniquify=False, useChirality=True, maxMatches=10001
    )
    if not 0 < len(matches) < 10001:
        raise ValueError("Incomplete independent stereochemical mapping")
    sdf_refs = bound.GetConformer().GetPositions()[np.array(matches)]
    rec = next(poses(d / "receptor.pdbqt"))
    tree = cKDTree([pos for typ, pos in rec.values() if typ not in ("HD", "HS", "H")])
    initial = next(poses(d / "ligand.pdbqt"))
    real = [i for i, (t, _) in initial.items() if t.startswith("CG")]
    dummies = [i for i, (t, _) in initial.items() if t.startswith("G")]
    if len(real) != 2 or len(dummies) != 2:
        raise ValueError("Expected a single closure edge")
    pairing = {
        i: min(dummies, key=lambda j: np.linalg.norm(initial[j][1] - initial[i][1]))
        for i in real
    }
    serial = meta["serial"]
    mapped = np.array(meta["mapped"])
    if (
        query.GetBondBetweenAtoms(
            next(i - 1 for i, s in order.items() if s == real[0]),
            next(i - 1 for i, s in order.items() if s == real[1]),
        )
        is None
    ):
        raise ValueError("Closure carbons do not share a chemical bond")
    deposited_distance = np.linalg.norm(
        mapped[:, serial.index(real[0]), :] - mapped[:, serial.index(real[1]), :],
        axis=1,
    )
    if not np.allclose(deposited_distance, deposited_distance[0]):
        raise ValueError("Ambiguous reference bond")
    return (
        serial,
        mapped,
        [order[i] for i in range(1, len(order) + 1)],
        sdf_refs,
        real,
        pairing,
        tree,
        deposited_distance[0],
    )


def main(work, output):
    with (output / "runs.csv").open() as handle:
        runs = list(csv.DictReader(handle))
    refs = {p: reference(work, p) for p in sorted({r["pdb"] for r in runs})}
    models = []
    for r in runs:
        serial, refs_a, order, refs_b, ends, pairs, tree, bond = refs[r["pdb"]]
        block = []
        for i, pose in enumerate(poses(output / "raw" / (r["key"] + ".pdbqt")), 1):
            coords = np.array([pose[s][1] for s in serial])
            xyz = np.array([pose[s][1] for s in order])
            rms_a = np.sqrt(((coords - refs_a) ** 2).sum(axis=-1).mean(axis=-1)).min()
            rms_b = np.sqrt(((xyz - refs_b) ** 2).sum(axis=-1).mean(axis=-1)).min()
            if abs(rms_a - rms_b) > 0.0015:
                raise ValueError(
                    f"{r['key']} {i}: deposited coordinate RMSDs differ: {rms_a} {rms_b}"
                )
            distances = tree.query(coords)[0]
            block.append(
                {
                    "key": r["key"],
                    "rank": i,
                    "rmsd_A": round(rms_a, 6),
                    "independent_sdf_rmsd_A": round(rms_b, 6),
                    "receptor_clash_atoms_2A": int((distances < 2).sum()),
                    "nearest_receptor_A": round(distances.min(), 6),
                    "cg_distance_A": round(
                        np.linalg.norm(pose[ends[0]][1] - pose[ends[1]][1]), 6
                    ),
                    "crystal_cg_distance_A": round(bond, 6),
                    "max_dummy_error_A": round(
                        max(
                            np.linalg.norm(pose[k][1] - pose[v][1])
                            for k, v in pairs.items()
                        ),
                        6,
                    ),
                }
            )
        if len(block) != int(r["modes"]) or not block:
            raise ValueError("Wrong mode count " + r["key"])
        top = block[0]
        scores = [m["rmsd_A"] for m in block]
        r.update(
            top1_A=top["rmsd_A"],
            best5_A=min(scores[:5]),
            oracle_A=min(scores),
            top1_clashes=top["receptor_clash_atoms_2A"],
            top1_cg_distance_A=top["cg_distance_A"],
            top1_endpoint_A=top["max_dummy_error_A"],
        )
        models.extend(block)
    write(output / "models.csv", models)
    write(output / "structures.csv", runs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    main(a.work, a.output)
