"""Compare the macrocyclic ring shape after a rigid fit; full-pose metrics stay in the receptor frame."""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np


def fit(moving, reference):
    a = moving - moving.mean(axis=0)
    b = reference - reference.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    u[:, -1] *= np.sign(np.linalg.det(u @ vt))
    return float(np.sqrt(np.mean(np.sum((a @ u @ vt - b) ** 2, axis=1))))


def main(a):
    sys.path.insert(0, str(a.repo.resolve() / "benchmarks/bond_distance"))
    from chem_graph import check

    sys.path.insert(0, str(a.repo.resolve() / "benchmarks/macrocycle_labels"))
    from geometry import coord_blocks
    from rdkit import Chem

    with (a.out / "runs.csv").open() as f:
        runs = list(csv.DictReader(f))
    rows = []
    for src in ("1FKJ", "6VCU"):
        m = json.loads((a.work / src / "preparation.json").read_text())
        inp_file = a.work / src / "ligand.pdbqt"
        graph, serial, ends, _, _, _ = check(
            inp_file, a.work / "raw" / (m["ligand"] + ".cif")
        )
        rings = Chem.GetSymmSSSR(graph)
        matches = [
            list(r)
            for r in rings
            if len(r) > 7
            and ends.issubset(
                {graph.GetAtomWithIdx(i).GetProp("source_name") for i in r}
            )
        ]
        if len(matches) != 1:
            raise ValueError("Ambiguous large ring spanning the cut")
        order = [{i: j for j, i in serial.items()}[i] for i in matches[0]]
        inp = next(coord_blocks(inp_file))
        initial = np.array([inp[j][1] for j in order])
        for r in runs:
            if r["query"] != src:
                continue
            receiver = r["receiver"]
            crystal = receiver if src == "1FKJ" else src
            ref = json.loads((a.work / crystal / "preparation.json").read_text())
            crystal_pos = [
                np.asarray(z)[[ref["serial"].index(i) for i in order]]
                for z in ref["mapped"]
            ]
            for rank, pose in enumerate(
                coord_blocks(a.out / "raw" / (r["key"] + ".pdbqt")), 1
            ):
                xyz = np.array([pose[i][1] for i in order])
                rows.append(
                    dict(
                        key=r["key"],
                        arm=r["arm"],
                        receiver=receiver,
                        query=src,
                        rank=rank,
                        ring_atoms=len(order),
                        input_to_crystal_A=min(fit(initial, z) for z in crystal_pos),
                        pose_to_input_A=fit(xyz, initial),
                        pose_to_crystal_A=min(fit(xyz, z) for z in crystal_pos),
                    )
                )
    with (a.out / "ring_shape.csv").open("w") as f:
        w = csv.DictWriter(f, rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("repo", "work", "out"):
        p.add_argument(k, type=Path)
    main(p.parse_args())
