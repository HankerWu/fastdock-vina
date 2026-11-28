"""Independently rebuild the BRD4 reference geometry and rerun every published-pose comparison."""

import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "cumulative/external"))
from audit import independent_pose, reference


def check(fixtures, results):
    cases = json.loads((fixtures / "manifest.json").read_text())
    for case in cases:
        for name, digest in case["sha256"].items():
            assert (
                hashlib.sha256((fixtures / case["pdb"] / name).read_bytes()).hexdigest()
                == digest
            )
    structures = {}
    for case in cases:
        if case["status"] != "ready":
            continue
        folder = fixtures / case["pdb"]
        info = json.loads((folder / "source.json").read_text())
        # The independent parser reads the first ligand occurrence and blank/A alternatives.
        lines = [
            l
            for l in (folder / (case["pdb"] + ".pdb")).read_text().splitlines()
            if l.startswith("HETATM")
            and l[17:20].strip() == case["ligand"]
            and l[21] == info["chain"]
            and l[76:78].strip() not in ("H", "D")
        ]
        assert lines[0][21:27] == info["ligand_instance"] and all(
            row[16] in " A" for row in lines if row[21:27] == info["ligand_instance"]
        )
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            (tmp / "raw").symlink_to(folder.resolve(), target_is_directory=True)
            structures[case["pdb"]] = reference(
                tmp, case["pdb"], case["ligand"], info["chain"]
            )
    ranks = list(csv.DictReader((results / "ranks.csv").open()))
    records = {(r["key"], int(r["rank"])): r for r in ranks}
    assert len(records) == len(ranks)
    maxerr, mappings, pocket, checked = 0.0, {}, {}, 0
    cells = list(csv.DictReader((results / "cells.csv").open()))
    for cell in cells:
        if cell["status"] != "ok":
            assert int(cell["modes"]) == 0
            continue
        q, r = cell["query"], cell["receiver"]
        qca, qcenter, mol = structures[q]
        rca, rcenter, _ = structures[r]
        ref = mol.GetConformer().GetPositions()
        if q != r:
            keys = sorted(
                k
                for k in qca.keys() & rca.keys()
                if np.linalg.norm(qca[k] - qcenter) < 14
                and np.linalg.norm(rca[k] - rcenter) < 14
            )
            assert len(keys) >= 12
            x, y = [np.asarray([ca[k] for k in keys]) for ca in (qca, rca)]
            u, _, v = np.linalg.svd((x - x.mean(0)).T @ (y - y.mean(0)))
            axis = np.diag([1.0, 1.0, np.linalg.det(u @ v)])
            rot = u @ axis @ v
            ref = (ref - x.mean(0)) @ rot + y.mean(0)
            pocket[q + "-" + r] = dict(
                n=len(keys),
                fit_rmsd=float(
                    np.mean(
                        np.sum(((x - x.mean(0)) @ rot + y.mean(0) - y) ** 2, axis=1)
                    )
                    ** 0.5
                ),
            )
        path = results / "raw" / (cell["key"] + ".pdbqt")
        assert hashlib.sha256(path.read_bytes()).hexdigest() == cell["pool_sha256"]
        observed = []
        for rank, block in enumerate(path.read_text().split("ENDMDL")[:-1], 1):
            distance, score, nmaps = independent_pose(block, mol, ref)
            row = records[cell["key"], rank]
            assert score == float(row["score"])
            mappings.setdefault(q, nmaps)
            assert mappings[q] == nmaps
            maxerr = max(maxerr, abs(distance - float(row["rmsd"])))
            observed.append(distance)
        assert len(observed) == int(cell["modes"])
        metrics = [
            observed[0],
            min(observed[:5]),
            min(observed),
            observed[0] - min(observed),
        ]
        maxerr = max(
            maxerr,
            max(
                abs(v - float(cell[f]))
                for v, f in zip(metrics, ["first", "best5", "best", "regret"])
            ),
        )
        checked += len(observed)
    assert checked == len(ranks) and maxerr < 1e-9
    return dict(
        cells=len(cells),
        published_poses=checked,
        max_difference_A=maxerr,
        stereo_mappings=mappings,
        pocket_alignments=pocket,
    )


if __name__ == "__main__":
    print(json.dumps(check(Path(sys.argv[1]), Path(sys.argv[2])), indent=2))
