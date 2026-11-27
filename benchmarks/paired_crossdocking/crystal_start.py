"""Diagnostic only: relax crystal starts and crystal torsions transplanted to the query geometry."""

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from run import alignment, coords


def write_ligand(mol, path):
    (setup,) = MoleculePreparation().prepare(mol)
    text, ok, err = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(err)
    path.write_text(text)
    return [
        l
        for l in text.splitlines()
        if l.startswith(("REMARK INDEX MAP", "BRANCH", "TORSDOF"))
    ]


def structures(sm, measured, tree):
    heavy = Chem.MolFromSmiles(sm["canonical_smiles"])
    true = Chem.Conformer(heavy.GetNumAtoms())
    for i, pt in enumerate(measured):
        true.SetAtomPosition(i, pt)
    heavy.AddConformer(true)
    crystal = Chem.AddHs(heavy, addCoords=True)
    start = Chem.AddHs(Chem.MolFromSmiles(sm["canonical_smiles"]))
    params = AllChem.ETKDGv3()
    params.randomSeed = 43103
    if AllChem.EmbedMolecule(start, params) or AllChem.MMFFOptimizeMolecule(
        start, maxIters=1000
    ):
        raise ValueError("Reference preparation changed")
    conf = start.GetConformer()
    indices = {s: i for i, s in enumerate(sm["serial"])}
    for line in tree.splitlines():
        if not line.startswith("BRANCH "):
            continue
        a, b = (indices[int(v)] for v in line.split()[1:])
        ta, tb = (heavy.GetAtomWithIdx(v) for v in (a, b))
        adj_a = sorted(n.GetIdx() for n in ta.GetNeighbors() if n.GetIdx() != b)
        adj_b = sorted(n.GetIdx() for n in tb.GetNeighbors() if n.GetIdx() != a)
        if not adj_a or not adj_b:
            continue
        torsion = (adj_a[0], a, b, adj_b[0])
        angle = rdMolTransforms.GetDihedralDeg(true, *torsion)
        rdMolTransforms.SetDihedralDeg(conf, *torsion, angle)
    points = conf.GetPositions()
    x = points[: heavy.GetNumAtoms()]
    xc, yc = x.mean(0), measured.mean(0)
    u, _, vt = np.linalg.svd((x - xc).T @ (measured - yc))
    u[:, -1] *= np.sign(np.linalg.det(u @ vt))
    for i, pt in enumerate((points - xc) @ u @ vt + yc):
        conf.SetAtomPosition(i, pt)
    return {"crystal": crystal, "projected": start}


def read_pose(path, serial):
    xyz = {
        int(l[6:11]): coords(l)
        for l in path.read_text().splitlines()
        if l.startswith(("ATOM  ", "HETATM"))
    }
    return np.array([xyz[k] for k in serial])


def rms(xyz, refs):
    return float(np.sqrt(np.sum((refs - xyz) ** 2, axis=2).mean(axis=1)).min())


def run(work, binary, pairs=None, compare_mappings=False):
    meta = json.loads((work / "metadata.json").read_text())
    folder = work / "crystal_diagnostic"
    folder.mkdir(exist_ok=True)
    rows = []
    projection_rows = []
    for source, sm in meta.items():
        for receiver, rm in meta.items():
            if sm["target"] != rm["target"] or (
                pairs is not None and (source, receiver) not in pairs
            ):
                continue
            rot, t, _, _ = alignment(sm, rm)
            refs = np.array(sm["mapped"]) @ rot + t
            template = (work / source / "ligand.pdbqt").read_text()
            starts = structures(sm, refs[0], template)
            if compare_mappings:
                # Equivalent atoms (e.g. sulfonyl oxygens) can define different
                # dihedral reference vectors in the prepared query. Choose by
                # initial heavy-atom fit before performing any relaxation.
                best = float("inf")
                for index, ref in enumerate(refs):
                    probe = (
                        starts["projected"]
                        if index == 0
                        else structures(sm, ref, template)["projected"]
                    )
                    xyz = probe.GetConformer().GetPositions()[: len(sm["serial"])]
                    initial = rms(xyz, refs)
                    projection_rows.append(
                        {
                            "source": source,
                            "receiver": receiver,
                            "mapping": index,
                            "initial_rmsd": initial,
                        }
                    )
                    if initial < best:
                        best = initial
                        starts["best_projected"] = probe
            for arm, mol in starts.items():
                base = f"{source}-{receiver}-{arm}"
                ligand, output = (
                    folder / (base + "-input.pdbqt"),
                    folder / (base + ".pdbqt"),
                )
                topology = write_ligand(mol, ligand)
                if topology != [
                    l
                    for l in template.splitlines()
                    if l.startswith(("REMARK INDEX MAP", "BRANCH", "TORSDOF"))
                ]:
                    raise ValueError(f"Topology changed for {base}")
                args = [
                    str(binary),
                    "--receptor",
                    str(work / receiver / "receptor.pdbqt"),
                    "--ligand",
                    str(ligand),
                    "--local_only",
                    "--seed",
                    "43103",
                    "--cpu",
                    "1",
                    "--out",
                    str(output),
                ]
                for axis, v in zip("xyz", rm["center"]):
                    args.extend(["--center_" + axis, str(v), "--size_" + axis, "24"])
                p = subprocess.run(args, check=True, text=True, capture_output=True)
                (folder / (base + ".log")).write_text(p.stdout + p.stderr)
                inter = float(
                    re.search(
                        r"Final Intermolecular Energy\s*:\s*([-\.\d]+)", p.stdout
                    )[1]
                )
                intra = float(
                    re.search(
                        r"Final Total Internal Energy\s*:\s*([-\.\d]+)", p.stdout
                    )[1]
                )
                rows.append(
                    {
                        "source": source,
                        "receiver": receiver,
                        "arm": arm,
                        "initial_rmsd": rms(read_pose(ligand, sm["serial"]), refs),
                        "rmsd": rms(read_pose(output, sm["serial"]), refs),
                        "inter_plus_intra": inter + intra,
                    }
                )
    if compare_mappings:
        with (work / "projection_mappings.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=projection_rows[0])
            writer.writeheader()
            writer.writerows(projection_rows)
    with (work / "crystal_starts.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    args = p.parse_args()
    run(args.work, args.binary.resolve())
