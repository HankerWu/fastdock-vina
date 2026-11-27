"""Probe the saved MAPK14 assay and BIRB-796 crystal geometry in two separate receptors."""

import argparse
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH / "pair_gradient/pose"))
from audit import coordinates, reference
from paired import RADII
from sources import crystal_contacts

spec = importlib.util.spec_from_file_location(
    "assay_docking", BENCH / "accepted_retention/affinity/run.py"
)
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)
IDS = {"1A9U": ("SB2", "176155"), "1KV2": ("B96", "156422")}


def ligand_identity(work, pdb, ccd, cid):
    ca, center, crystal = reference(work, pdb, ccd, "A")
    query, _, _ = coordinates(
        (work / "prepared" / (cid + ".pdbqt")).read_text().splitlines()
    )
    shape = Chem.Mol(query)
    Chem.AssignStereochemistryFrom3D(shape, replaceExistingTags=True)
    matches = crystal.GetSubstructMatches(
        shape, uniquify=False, useChirality=True, maxMatches=10001
    )
    assert len(matches) < 10001
    identity = {
        "pdb": pdb,
        "component": ccd,
        "cid": cid,
        "heavy_atoms": query.GetNumAtoms(),
        "query_input_smiles": Chem.MolToSmiles(query),
        "query_3d_stereo": Chem.MolToSmiles(shape),
        "crystal_3d_stereo": Chem.MolToSmiles(crystal),
        "input_3d_inchikey": Chem.MolToInchiKey(shape),
        "crystal_inchikey": Chem.MolToInchiKey(crystal),
        "stereo_matches": len(matches),
        "topology_matches": len(
            crystal.GetSubstructMatches(query, uniquify=False, useChirality=False)
        ),
    }
    return ca, center, crystal, shape, matches, identity


def rcoords(line):
    return [float(line[i : i + 8]) for i in (30, 38, 46)]


def relax_crystal(work, pdb, parent, probe, projected, maps):
    from crystal_start import write_ligand

    folder = work / "local"
    folder.mkdir(exist_ok=True)
    ligand, output = folder / (pdb + "-input.pdbqt"), folder / (pdb + "-relaxed.pdbqt")
    mol = Chem.Mol(probe)
    mol.RemoveAllConformers()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, point in enumerate(projected[list(maps[0])]):
        conf.SetAtomPosition(i, point)
    mol.AddConformer(conf)
    write_ligand(Chem.AddHs(mol, addCoords=True), ligand)

    # Atom serials can be reordered within rigid fragments by Meeko. Compare
    # directed torsions and heavy atom types/charges in canonical SMILES order.
    def topology(path):
        lines = path.read_text().splitlines()
        idx = {}
        for line in lines:
            if line.startswith("REMARK SMILES IDX"):
                nums = list(map(int, line.split()[3:]))
                idx.update(zip(nums[1::2], nums[::2], strict=True))
        bonds = set()
        atoms = {}
        for line in lines:
            if line.startswith("BRANCH "):
                a, b = map(int, line.split()[1:])
                bonds.add((idx[a], idx[b]))
            if line.startswith(("ATOM  ", "HETATM")):
                serial = int(line[6:11])
                if serial in idx:
                    atoms[idx[serial]] = tuple(line.split()[-2:])
        return bonds, atoms

    if topology(ligand) != topology(work / "prepared/156422.pdbqt"):
        raise ValueError(
            "Crystal and query graphs disagree on torsions or atom parameters"
        )
    command = [
        parent,
        "--receptor",
        work / "prepared" / (pdb + ".pdbqt"),
        "--ligand",
        ligand,
        "--config",
        work / "prepared" / (pdb + ".box"),
        "--cpu",
        "1",
        "--seed",
        "458921",
        "--local_only",
        "--out",
        output,
    ]
    run = subprocess.run(command, text=True, capture_output=True, check=True)
    (folder / (pdb + ".log")).write_text(run.stdout + run.stderr)
    _, xyz, _ = coordinates(output.read_text().splitlines())
    return min(
        float(np.sqrt(np.square(xyz - projected[list(m)]).sum(1).mean())) for m in maps
    )


def analyze(work, parent, quadratic, halving):
    meta = {p: ligand_identity(work, p, c, cid) for p, (c, cid) in IDS.items()}
    (work / "identity.json").write_text(
        json.dumps([m[-1] for m in meta.values()], indent=2) + "\n"
    )
    qca, qcenter, crystal, probe, maps, _ = meta["1KV2"]
    if not maps or probe.GetNumAtoms() != crystal.GetNumAtoms():
        raise ValueError("BIRB-796 must match the full query graph")
    if meta["1A9U"][4]:
        raise ValueError(
            "SB-203580 stereo status differs from the frozen input; reassess the eligibility"
        )
    qx = crystal.GetConformer().GetPositions()
    radii = np.array([RADII[a.GetSymbol()] for a in crystal.GetAtoms()])
    pre = work / "prepared"
    pooled = engine.read(work / "observed/runs.csv")
    modes = []
    cells = []
    geometries = []
    for pdb, (component, _) in IDS.items():
        ca, center = meta[pdb][:2]
        keys = sorted(
            k
            for k in ca.keys() & qca.keys()
            if np.linalg.norm(qca[k] - qcenter) < 14
            and np.linalg.norm(ca[k] - center) < 14
        )
        q, r = np.array([qca[k] for k in keys]), np.array([ca[k] for k in keys])
        rotation, _ = Rotation.align_vectors(r - r.mean(0), q - q.mean(0))
        projected = rotation.apply(qx - q.mean(0)) + r.mean(0)
        residuals = rotation.apply(q - q.mean(0)) + r.mean(0) - r
        rec = (pre / (pdb + ".pdbqt")).read_text().splitlines()
        heavy = [
            s
            for s in rec
            if s.startswith(("ATOM", "HETATM")) and s.split()[-1] not in ("HD", "H")
        ]
        d = cdist(projected, np.array([rcoords(s) for s in heavy]))
        limit = 0.75 * (
            radii[:, None] + np.array([RADII[s.split()[-1]] for s in heavy])
        )
        clashes = np.sum(d < limit, axis=0)
        from collections import Counter

        residues = Counter()
        for line, n in zip(heavy, clashes, strict=True):
            if n:
                residues[line[17:26].strip()] += int(n)
        box = np.array(
            [
                float(l.split("=")[1])
                for l in (pre / (pdb + ".box")).read_text().splitlines()
                if l.startswith("center_")
            ]
        )
        geometry = {
            "pdb": pdb,
            "pocket_CA": len(keys),
            "pocket_CA_RMSD_A": float(np.sqrt((residuals**2).sum(1).mean())),
            "closest_A": float(d.min()),
            "reference_clashes": int(clashes.sum()),
            "clashes_by_residue": dict(residues),
            "out_of_box_atoms": int((abs(projected - box) > 13).any(1).sum()),
            "center_shift_A": float(np.linalg.norm(projected.mean(0) - box)),
        }
        raw = (work / "raw" / (pdb + ".pdb")).read_text().splitlines()
        crystal_struct = gemmi.read_structure(str(work / "raw" / (pdb + ".pdb")))
        geometry.update(crystal_contacts(crystal_struct, "A", component))
        dep = [
            s
            for s in raw
            if s.startswith("ATOM")
            and s[21] == "A"
            and s[16] in " A"
            and s[76:78].strip() not in ("H", "D")
        ]
        original = {(s[12:16].strip(), s[22:27].strip()) for s in dep}
        rebuilt = [
            s for s in heavy if (s[12:16].strip(), s[22:27].strip()) not in original
        ]
        geometry["new_heavy_atoms_within_8A"] = [
            s[12:27]
            for s in rebuilt
            if np.linalg.norm(projected - rcoords(s), axis=1).min() < 8
        ]
        geometry["min_repaired_to_BIRB_A"] = min(
            (
                float(np.linalg.norm(projected - rcoords(s), axis=1).min())
                for s in rebuilt
            ),
            default=None,
        )
        geometry["resolved_gap_flanks_A"] = {
            str(i): float(
                min(
                    np.linalg.norm(projected - rcoords(s), axis=1).min()
                    for s in dep
                    if s[22:26].strip() == str(i)
                )
            )
            for i in (114, 123, 169, 185)
        }
        geometry["alignment_sensitivity"] = {}
        for label, common in (
            (
                "pocket_without_activation_loop",
                [k for k in keys if not 165 <= k[0] <= 190],
            ),
            ("whole_chain", sorted(ca.keys() & qca.keys())),
        ):
            x, y = np.array([qca[k] for k in common]), np.array([ca[k] for k in common])
            rot, _ = Rotation.align_vectors(y - y.mean(0), x - x.mean(0))
            transformed = rot.apply(qx - x.mean(0)) + y.mean(0)
            dist = cdist(transformed, [rcoords(line) for line in heavy])
            phe = [
                j
                for j, line in enumerate(heavy)
                if line[17:20] == "PHE" and line[22:26].strip() == "169"
            ]
            geometry["alignment_sensitivity"][label] = {
                "n_CA": len(common),
                "clashes": int((dist < limit).sum()),
                "Phe169_clashes": int((dist[:, phe] < limit[:, phe]).sum()),
                "closest_A": float(dist.min()),
            }
        geometry["crystal_start_local_rmsd_A"] = relax_crystal(
            work, pdb, parent, probe, projected, maps
        )
        geometry["resolution_A"] = crystal_struct.resolution
        geometries.append(geometry)
        for i, seed in enumerate((458921, 703109)):
            arms = ("parent", "quadratic") if i == 0 else ("quadratic", "parent")
            for arm in (*arms, "halving"):
                if pdb == "1A9U":
                    row = next(
                        r for r in pooled if int(r["seed"]) == seed and r["arm"] == arm
                    )
                    folder = work / "observed"
                else:
                    cell = (
                        "p38-alpha",
                        pdb,
                        "156422",
                        pre / (pdb + ".pdbqt"),
                        pre / (pdb + ".box"),
                        pre / "156422.pdbqt",
                    )
                    row = engine.dock(
                        cell,
                        work / "run",
                        {"parent": parent, "quadratic": quadratic, "halving": halving}[
                            arm
                        ],
                        arm,
                        seed,
                        12000,
                        True,
                    )
                    folder = work / "run/raw"
                vals = []
                for j, block in enumerate(
                    (folder / (row["key"] + ".pdbqt")).read_text().split("ENDMDL")[:-1],
                    1,
                ):
                    mol, xyz, _ = coordinates(block.splitlines())
                    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(probe)
                    rms = min(
                        float(
                            np.sqrt(
                                np.square(xyz - projected[np.array(m)]).sum(1).mean()
                            )
                        )
                        for m in maps
                    )
                    vals.append(rms)
                    dist = cdist(xyz, [rcoords(line) for line in heavy])
                    threshold = limit[np.array(maps[0])]
                    modes.append(
                        {
                            "key": row["key"],
                            "rank": j,
                            "rmsd_A": rms,
                            "clashes": int((dist < threshold).sum()),
                            "closest_A": float(dist.min()),
                        }
                    )
                cells.append(
                    {
                        **row,
                        "top1_A": vals[0],
                        "top5_A": min(vals[:5]),
                        "oracle_A": min(vals),
                    }
                )
    (work / "geometry.json").write_text(json.dumps(geometries, indent=2) + "\n")
    engine.write(work / "cells.csv", cells)
    engine.write(work / "modes.csv", modes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("parent", type=Path)
    parser.add_argument("quadratic", type=Path)
    parser.add_argument("halving", type=Path)
    parser.add_argument("--cpu", type=int, required=True)
    a = parser.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    (a.work / "run/raw").mkdir(parents=True, exist_ok=True)
    analyze(
        a.work.resolve(), a.parent.resolve(), a.quadratic.resolve(), a.halving.resolve()
    )
