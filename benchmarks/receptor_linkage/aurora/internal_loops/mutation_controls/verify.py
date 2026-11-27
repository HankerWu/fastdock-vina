"""Independent receptor-frame, bond, stereo and empirical-quantile checks."""

import argparse
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from rdkit import Chem

from run import inputs, measure, positions, replay


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def audit(aurora, brd, controls, scan, out):
    out.mkdir(exist_ok=False)
    with TemporaryDirectory(dir=out) as temp:
        systems = inputs(aurora, brd, Path(temp))
        producer = json.loads((scan / "inputs.json").read_text())
        assert all(
            all(
                replay.digest(p) == producer[name][k]
                for p, k in zip(item[:3], ["receptor", "ligand", "bonds"], strict=True)
            )
            for name, item in systems.items()
        )
    errors = []
    bonds = []
    cells = []
    n_stereo = 0
    n_modes = 0
    for result, key, folder in ((controls, "ligand", "pools"), (scan, "drug", "raw")):
        reported = {
            (r[key], r["arm"], r["seed"], int(r["rank"])): r
            for r in rows(result / "modes.csv")
        }
        cell_modes = 0
        for c in rows(result / "cells.csv"):
            name, arm, seed = c[key], c["arm"], c["seed"]
            path = result / folder / f"{name}_{seed}_{arm}.pdbqt"
            assert replay.digest(path) == c["pool_sha256"]
            _, ligand, _, _, _, mol, serials, ref = systems[name]
            measured = measure(path, ligand, mol, serials, ref)
            assert len(measured) == int(c["n"])
            cell_modes += len(measured)
            for rank, (r, score, drift) in enumerate(measured, 1):
                row = reported[name, arm, seed, rank]
                errors.append(abs(r - float(row["rmsd_A"])))
                bonds.append(drift)
                assert score == float(row["score"])
            # Compare signed volumes at specified stereocentres, in input SMILES
            # order. For this panel only danusertib has a specified centre (CH).
            lines = ligand.read_text().splitlines()
            query = Chem.MolFromSmiles(
                next(i.split()[2] for i in lines if i.startswith("REMARK SMILES "))
            )
            centers = Chem.FindMolChiralCenters(query, includeUnassigned=True)
            template = {
                int(i[6:11]): [float(i[k : k + 8]) for k in (30, 38, 46)]
                for i in lines
                if i.startswith(("ATOM  ", "HETATM"))
            }
            for center, tag in centers:
                assert tag != "?"
                n = [
                    serials[a.GetIdx()]
                    for a in query.GetAtomWithIdx(center).GetNeighbors()
                ]
                assert len(n) == 3

                def sign(x):
                    return np.sign(
                        np.linalg.det(
                            [np.subtract(x[i], x[serials[center]]) for i in n]
                        )
                    )

                assert sign(template) != 0
                assert all(sign(p) == sign(template) for p in positions(path))
                n_stereo += len(measured)
            vals = [v[0] for v in measured]
            summary = dict(
                top1=vals[0],
                best5=min(vals[:5]),
                oracle=min(vals),
                regret=vals[0] - min(vals),
            )
            errors.extend(abs(v - float(c[k])) for k, v in summary.items())
            cells.append(
                dict(
                    ligand=name,
                    seed=int(seed),
                    arm=arm,
                    **summary,
                    n=len(vals),
                    cpu_s=float(c["cpu_s"]),
                    wall_s=float(c["wall_s"]),
                )
            )
        assert cell_modes == len(reported)
        n_modes += cell_modes
    assert max(errors) < 1e-8
    summary = []
    for arm in ("native", "scan", "official", "pro_double"):
        for drug in ("all", "at9283", "danusertib", "82I"):
            selected = [
                r
                for r in cells
                if r["seed"] != 61231
                and r["arm"] == arm
                and drug in ("all", r["ligand"])
            ]
            assert len(selected) == (6 if drug == "all" else 2)
            top = np.asarray([r["top1"] for r in selected])
            summary.append(
                dict(
                    arm=arm,
                    ligand=drug,
                    n=len(top),
                    mean=top.mean(),
                    **dict(
                        zip(
                            ("q25", "median", "q75", "p90", "p95"),
                            np.percentile(top, [25, 50, 75, 90, 95]),
                        )
                    ),
                    **{
                        k: np.mean([r[k] for r in selected])
                        for k in ("best5", "oracle", "regret", "cpu_s", "wall_s")
                    },
                    severe_ge8=sum(top >= 8),
                    empty=sum(r["n"] == 0 for r in selected),
                )
            )
    replay.geom.write_csv(out, "distributions.csv", summary)
    summary = dict(
        cells=len(cells),
        modes=n_modes,
        maximum_RMSD_difference_A=max(errors),
        maximum_heavy_bond_drift_A=max(bonds),
        checked_signed_volumes=n_stereo,
        sign_changes=0,
        empty_pools=0,
        matched_input_triplets=len(systems),
    )
    (out / "verification.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("aurora", "brd", "controls", "scan", "out"):
        p.add_argument(k, type=Path)
    audit(*(v.resolve() for v in vars(p.parse_args()).values()))
