"""Paired donor starts on alternate crystal structures; test labels enter only after docking."""

import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from crystal_start import read_pose, rms, write_ligand
from rdkit import Chem
from rdkit.Chem import AllChem
from run import ROOT, alignment

sys.path.insert(0, str(ROOT / "tools"))
from make_donor_starts import make_starts


def table(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def run(work, binary):
    folder = work / "donor_starts"
    folder.mkdir(exist_ok=True)
    meta = json.loads((work / "metadata.json").read_text())
    rows = []
    cells = []
    for s, sm in meta.items():
        for r, rm in meta.items():
            if s == r or sm["target"] != rm["target"]:
                continue
            q = Chem.AddHs(Chem.MolFromSmiles(sm["canonical_smiles"]))
            ep = AllChem.ETKDGv3()
            ep.randomSeed = 43103
            if AllChem.EmbedMolecule(q, ep) or AllChem.MMFFOptimizeMolecule(
                q, maxIters=1000
            ):
                raise ValueError("query embedding")
            donor = Chem.MolFromSmiles(rm["canonical_smiles"])
            dc = Chem.Conformer(donor.GetNumAtoms())
            for i, pt in enumerate(rm["mapped"][0]):
                dc.SetAtomPosition(i, pt)
            donor.AddConformer(dc)
            for copy in (False, True):
                arm = "transfer" if copy else "align"
                candidates = make_starts(q, donor, copy)
                for k, (core_error, mapping, mol) in enumerate(candidates):
                    base = f"{s}-{r}-{arm}-{k}"
                    inp, out = (
                        folder / (base + "-input.pdbqt"),
                        folder / (base + ".pdbqt"),
                    )
                    tree = write_ligand(mol, inp)
                    if tree != [
                        l
                        for l in (work / s / "ligand.pdbqt").read_text().splitlines()
                        if l.startswith(("REMARK INDEX MAP", "BRANCH", "TORSDOF"))
                    ]:
                        raise ValueError("Changed PDBQT topology")
                    cmd = [
                        str(binary),
                        "--receptor",
                        str(work / r / "receptor.pdbqt"),
                        "--ligand",
                        str(inp),
                        "--out",
                        str(out),
                        "--local_only",
                        "--cpu",
                        "1",
                        "--seed",
                        "43103",
                    ]
                    for axis, c in zip("xyz", rm["center"]):
                        cmd.extend(["--center_" + axis, str(c), "--size_" + axis, "24"])
                    log = subprocess.run(
                        cmd, check=True, text=True, capture_output=True
                    ).stdout
                    (folder / (base + ".log")).write_text(log)
                    score = sum(
                        float(re.search(pattern + r"\s*:\s*([-\.\d]+)", log)[1])
                        for pattern in (
                            "Final Intermolecular Energy",
                            "Final Total Internal Energy",
                        )
                    )
                    # Query crystal geometry is used only here, after the output exists.
                    rot, t, _, _ = alignment(sm, rm)
                    refs = np.array(sm["mapped"]) @ rot + t
                    rows.append(
                        {
                            "target": sm["target"],
                            "source": s,
                            "receiver": r,
                            "arm": arm,
                            "candidate": k,
                            "mapped_atoms": len(mapping),
                            "query_heavy": q.GetNumHeavyAtoms(),
                            "core_rmsd": core_error,
                            "inter_plus_intra": score,
                            "rmsd": rms(read_pose(out, sm["serial"]), refs),
                        }
                    )
    for path in sorted(work.glob("poses_*.csv")):
        cap = path.stem.split("_")[1]
        for row in csv.DictReader(path.open()):
            if row["kind"] != "cross" or row["rank"] != "1" or row["arm"] != "native":
                continue
            s, r = row["ligand"], row["receiver"]
            pool = work / ("runs_" + cap) / (f"{s}_{r}_{row['seed']}_native.pdbqt")
            e = float(
                re.search(r"REMARK INTER \+ INTRA:\s*([-\.\d]+)", pool.read_text())[1]
            )
            for arm in ["align", "transfer"]:
                extra = [
                    x
                    for x in rows
                    if (x["source"], x["receiver"], x["arm"]) == (s, r, arm)
                ]
                scale = float(row["rmsd"])
                for x in sorted(extra, key=lambda x: x["inter_plus_intra"])[:1]:
                    if x["inter_plus_intra"] < e:
                        scale = x["rmsd"]
                cells.append(
                    {
                        "target": row["target"],
                        "source": s,
                        "receiver": r,
                        "cap": cap,
                        "seed": row["seed"],
                        "arm": arm,
                        "n_added": len(extra),
                        "native_rmsd": row["rmsd"],
                        "selected_rmsd": scale,
                        "native_energy": e,
                    }
                )
    table(work / "donor_poses.csv", rows)
    table(work / "donor_cells.csv", cells)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    a = p.parse_args()
    run(a.work, a.binary.resolve())
