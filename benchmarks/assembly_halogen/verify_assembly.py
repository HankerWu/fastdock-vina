"""Independently verify assembly pose endpoints from BIOMT and Meeko conformers."""

import argparse
import csv
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate
from verify import graph_maps


def computed_site_refs(folder, mol):
    ref = Chem.SDMolSupplier(str(folder / "crystal.sdf"), removeHs=True)[0]
    x = ref.GetConformer().GetPositions()
    lines = (folder / "deposition.pdb").read_text().splitlines()
    # All follow-up cases have a single biological unit and transformations that act on the ligand chain.
    operations = {}
    for line in lines:
        if line.startswith("REMARK 350   BIOMT"):
            fields = line.split()
            k = int(fields[2][-1]) - 1
            index = int(fields[3])
            operations.setdefault(index, np.zeros((3, 4)))[k] = [
                float(v) for v in fields[4:8]
            ]
    coordinates = [x]
    for transform in operations.values():
        xyz = x @ transform[:, :3].T + transform[:, 3]
        if 1e-7 < np.linalg.norm(xyz.mean(0) - x.mean(0)) <= 2:
            coordinates.append(xyz)
    maps = np.array(graph_maps(mol, ref, False))
    return [xyz[maps] for xyz in coordinates]


def verify(work, assembly, results):
    rows = []
    for name in ("assembly-modes.csv", "assembly-modes-high.csv"):
        with (results / name).open() as f:
            expected = list(csv.DictReader(f))
        paths = {}
        for r in expected:
            case = r["case"]
            seed = int(r["seed"])
            arm = r["arm"]
            cap = int(r["cap"])
            path = (
                assembly / case / f"{arm}-{seed}-cap{cap}.pdbqt"
                if "high" in name
                else assembly / case / f"{arm}-{seed}.pdbqt"
            )
            if not path.exists() and arm == "asu":
                path = work / case / f"poses-{seed}.pdbqt"
            paths.setdefault(path, []).append(r)
        max_error = 0
        models = 0
        for path, targets in paths.items():
            m = Chem.RemoveHs(
                RDKitMolCreate.from_pdbqt_mol(
                    PDBQTMolecule.from_file(str(path), skip_typing=True)
                )[0]
            )
            refs = computed_site_refs(work / targets[0]["case"], m)
            for c, r in zip(m.GetConformers(), targets, strict=True):
                d = [
                    np.sqrt(
                        np.square(array - c.GetPositions()).sum(axis=2).mean(axis=1)
                    ).min()
                    for array in refs
                ]
                max_error = max(
                    max_error,
                    abs(float(r["rmsd_A"]) - d[0]),
                    abs(float(r["site_rmsd_A"]) - min(d)),
                )
                models += 1
        rows.append(
            dict(
                source=name, pools=len(paths), models=models, max_RMSD_error_A=max_error
            )
        )
        # BIOMTs operate on rounded input coordinates; official coordinate files round again.
        if max_error > 0.002:
            raise ValueError(f"{name}: coordinates differ {max_error}")
    (results / "assembly-independent.json").write_text(
        json.dumps(rows, indent=2) + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("assembly", type=Path)
    p.add_argument("results", type=Path)
    a = p.parse_args()
    verify(a.work, a.assembly, a.results)
