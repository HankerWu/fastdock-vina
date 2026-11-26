"""Reconstruct receptor-frame RMSDs from independent ModelServer bound SDFs, raw protein and Meeko."""

import argparse
import csv
import hashlib
import json
import urllib.request
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

HERE = Path(__file__).resolve().parent


def fit_input(work, p):
    meta = json.loads((work / p / "source.json").read_text())
    rows = (work / p / (p + ".pdb")).read_text().splitlines()
    ident = meta["ligand_instance"]
    url = f"https://models.rcsb.org/v1/{p.lower()}/ligand?auth_asym_id={ident[0]}&auth_seq_id={int(ident[1:5])}&encoding=sdf"
    dest = work / p / "observed.sdf"
    if not dest.exists():
        with urllib.request.urlopen(url, timeout=60) as response:
            dest.write_bytes(response.read())
    ligand = Chem.SDMolSupplier(str(dest), removeHs=True)[0]
    xyz = ligand.GetConformer().GetPositions()
    ca = {}
    for line in rows:
        if (
            line.startswith("ATOM  ")
            and line[21] == meta["chain"]
            and line[12:16].strip() == "CA"
            and line[16] in (" ", "A")
        ):
            ca[(line[22:27], line[17:20])] = [
                float(line[k : k + 8]) for k in (30, 38, 46)
            ]
    return (
        ligand,
        xyz.mean(axis=0),
        ca,
        {
            k: v
            for k, v in [
                ("url", url),
                ("sha256", hashlib.sha256(dest.read_bytes()).hexdigest()),
            ]
        },
    )


def verify(work, out):
    expected = list(csv.DictReader((out / "modes.csv").open()))
    cases = {}
    sources = {}
    for pdb in sorted({r["query"] for r in expected}):
        lig, center, ca, src = fit_input(work, pdb)
        cases[pdb] = (lig, center, ca)
        sources[pdb] = src
    checks = []
    for file in sorted((work / "pools").glob("*.pdbqt")):
        q, r, seed = file.stem.split("-")
        qlig, xmid, x = cases[q]
        _, ymid, y = cases[r]
        keys = sorted(
            k
            for k in x.keys() & y.keys()
            if np.linalg.norm(x[k] - xmid) < 14 and np.linalg.norm(y[k] - ymid) < 14
        )
        xx, yy = [np.asarray([c[k] for k in keys]) for c in (x, y)]
        transform, _ = Rotation.align_vectors(yy - yy.mean(0), xx - xx.mean(0))
        ref = Chem.Mol(qlig)
        coords = transform.apply(
            qlig.GetConformer().GetPositions() - xx.mean(0)
        ) + yy.mean(0)
        for i, one in enumerate(coords):
            ref.GetConformer().SetAtomPosition(i, one)
        (m,) = RDKitMolCreate.from_pdbqt_mol(
            PDBQTMolecule.from_file(str(file), skip_typing=True)
        )
        m = Chem.RemoveHs(m)
        assert Chem.MolToSmiles(m) == Chem.MolToSmiles(ref), (
            q,
            Chem.MolToSmiles(m),
            Chem.MolToSmiles(ref),
        )
        base = [
            v
            for v in expected
            if (v["query"], v["receiver"], v["seed"]) == (q, r, seed)
        ]
        assert len(base) == m.GetNumConformers()
        for i, v in enumerate(base):
            dist = rdMolAlign.CalcRMS(m, ref, prbId=i, refId=0)
            checks.append(
                dict(
                    query=q,
                    receiver=r,
                    seed=seed,
                    model=i + 1,
                    ca_count=len(keys),
                    rmsd_A=dist,
                    abs_difference_A=abs(dist - float(v["rmsd_A"])),
                )
            )
    assert len(checks) == len(expected)
    err = max(c["abs_difference_A"] for c in checks)
    if err > 0.001:
        raise ValueError(f"Independent RMSD discrepancy {err}")
    (out / "sdf_sources.json").write_text(json.dumps(sources, indent=2) + "\n")
    with (out / "sdf_check.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=checks[0])
        w.writeheader()
        w.writerows(checks)
    print(len(checks), "independent SDF RMSDs, max absolute error", err)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--out", type=Path, default=HERE)
    a = p.parse_args()
    verify(a.work, a.out)
