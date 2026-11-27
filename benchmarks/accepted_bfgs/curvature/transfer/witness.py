"""Post hoc feasibility probe. Fit legal PDBQT torsions to the transferred crystal."""

import argparse
import itertools
import re
import json
import sys
import subprocess
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from evaluate import PDBS, frame, write_csv

BENCH = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(BENCH / "receptor_linkage/aurora"))
from pose_probe import (
    get_tree,
    fit_torsions,
    local_refinement,
    mapped_rmsd,
    register,
    score_pose,
    store_pose,
    run_logged,
)
from geometry import mode_coords


def probe(work, parent, fixed, parser):
    meta = {p: json.loads((work / p / "metadata.json").read_text()) for p in PDBS}
    out = work / "witness"
    out.mkdir(exist_ok=True)
    rows, amides = [], []
    for q, r in itertools.permutations(PDBS):
        path = work / q / "ligand.pdbqt"
        ref, _, _ = frame(meta[q], meta[r])
        mol = Chem.MolFromSmiles(meta[q]["smiles"])
        complete = Chem.AddHs(mol)
        lines, ids, order, coords, arms = get_tree(path)
        index = {v: i for i, v in enumerate(order)}
        heavy = np.array([index[i] for i in range(mol.GetNumAtoms())])
        hset = set(heavy)
        visible = [arm for arm in arms if len(hset.intersection(arm[2]) - {arm[1]})]
        # Construct one admissible witness. Using the first graph isomorphism
        # avoids treating unconverged local fits as global minima. Final
        # geometric evaluation still checks every stereochemical isomorphism.
        best, fit, _ = fit_torsions(
            coords, heavy, visible, ref[:1], accept_finite_limit=True
        )
        placed = out / (q + "_fit.pdbqt")
        store_pose(placed, lines, fit)
        edges = {
            tuple(
                sorted((ids[index[b.GetBeginAtomIdx()]], ids[index[b.GetEndAtomIdx()]]))
            )
            for b in complete.GetBonds()
            if b.GetBeginAtomIdx() in index and b.GetEndAtomIdx() in index
        }
        bonds = out / (q + ".bonds")
        bonds.write_text("".join(f"{a} {b}\n" for a, b in sorted(edges)))
        signature = subprocess.check_output([parser, path, bonds])
        inferred = subprocess.check_output([parser, path])
        matches_native = inferred == signature
        (out / f"{q}_explicit_topology.txt").write_bytes(signature)
        (out / f"{q}_prepared_topology.txt").write_bytes(inferred)
        assert subprocess.check_output([parser, placed, bonds]) == signature
        positions = {sid: i for i, sid in enumerate(ids)}
        invariant_pairs = {(positions[x], positions[y]) for x, y in edges}
        for atom in complete.GetAtoms():
            if atom.GetIdx() not in index:
                continue
            neighbors = [
                index[n.GetIdx()] for n in atom.GetNeighbors() if n.GetIdx() in index
            ]
            invariant_pairs.update(itertools.combinations(neighbors, 2))

        def bond_change(moved):
            return max(
                abs(
                    np.linalg.norm(coords[a] - coords[b])
                    - np.linalg.norm(moved[a] - moved[b])
                )
                for a, b in invariant_pairs
            )

        assert bond_change(fit) < 1e-10
        amide_bonds = {
            tuple(sorted(indices[1:3]))
            for indices in mol.GetSubstructMatches(Chem.MolFromSmarts("[O]=[C]-[N]"))
        }
        restricted = [
            arm
            for arm in visible
            if tuple(sorted((order[arm[0]], order[arm[1]]))) not in amide_bonds
        ]
        _, fit_locked, _ = fit_torsions(
            coords, heavy, restricted, ref[:1], accept_finite_limit=True
        )
        assert bond_change(fit_locked) < 1e-10
        locked_path = out / f"{q}_all_amides_locked.pdbqt"
        store_pose(locked_path, lines, fit_locked)
        rotbonds = {tuple(sorted((int(order[a]), int(order[b])))) for a, b, _ in arms}
        sel = int(np.argmin(np.mean(np.sum((ref - fit[heavy]) ** 2, axis=2), axis=1)))
        tmpl = ref[sel]
        smarts = Chem.MolFromSmarts("[O]=[C]-[N]-[C]")
        for indices in mol.GetSubstructMatches(smarts):
            values = []
            for xyz in (coords[heavy], fit[heavy], tmpl):
                cf = Chem.Conformer(len(heavy))
                for i, v in enumerate(xyz):
                    cf.SetAtomPosition(i, v)
                values.append(rdMolTransforms.GetDihedralDeg(cf, *indices))
            amides.append(
                dict(
                    query=q,
                    indices_zero_based=str(indices),
                    axis_mobile=tuple(sorted(indices[1:3])) in rotbonds,
                    prepared_deg=values[0],
                    fitted_deg=values[1],
                    mapped_crystal_deg=values[2],
                )
            )
        common = [
            "--receptor",
            work / r / "receptor.pdbqt",
            "--seed",
            60913,
            "--cpu",
            1,
        ]
        for axis, value in zip("xyz", meta[r]["center"]):
            common += ["--center_" + axis, value, "--size_" + axis, 24]
        for arm, exe in [("parent", parent), ("fixed", fixed)]:
            localized = out / f"{q}_{arm}_local.pdbqt"
            local_refinement(
                [exe, *common], placed, bonds, localized, localized.with_suffix(".log")
            )
            scoring = [exe, *common, "--ligand", path, "--ligand_bonds", bonds]
            stages = {
                "fit": placed,
                "all_amides_locked": locked_path,
                "local": localized,
            }
            for se in [27183, 91439]:
                winner = out / f"{q}_{arm}_{se}.pdbqt"
                first = (
                    (work / "runs" / f"{q}_{r}_{se}_{arm}.pdbqt")
                    .read_text()
                    .split("ENDMDL", 1)[0]
                )
                winner.write_text(first + "ENDMDL\n")
                stages[str(se)] = winner
            for stage, f in stages.items():
                e = score_pose(scoring, f)
                original = run_logged(
                    [
                        exe,
                        *common,
                        "--ligand",
                        path,
                        "--pose_coordinates",
                        f,
                        "--score_only",
                        "--unbound_energy",
                        0,
                    ],
                    out / f"{q}_{arm}_{stage}_inferred_score.log",
                )
                implicit_e = float(
                    re.search(
                        r"Estimated Free Energy of Binding\s*:\s*([-+.\d]+)", original
                    )[1]
                )
                rows.append(
                    dict(
                        query=q,
                        receiver=r,
                        arm=arm,
                        stage=stage,
                        score_unbound_zero=e,
                        inferred_score_unbound_zero=implicit_e,
                        rmsd=mapped_rmsd(mode_coords(f, meta[q]["serial"]), ref),
                        rigid_only=min(
                            np.linalg.norm(register(coords, heavy, t)[heavy] - t)
                            / np.sqrt(len(heavy))
                            for t in ref
                        ),
                        fit_first_mapping_rmsd=best,
                        all_branches=len(arms),
                        heavy_visible=len(visible),
                        n_ccd_bonds_retained=len(edges),
                        inferred_equals_ccd=matches_native,
                        constructed_rmsd_before_rounding=mapped_rmsd(fit[heavy], ref),
                        locked_rmsd_before_rounding=mapped_rmsd(fit_locked[heavy], ref),
                        fitted_max_bond_or_one_three_change=bond_change(fit),
                        amide_locked_max_bond_or_one_three_change=bond_change(
                            fit_locked
                        ),
                    )
                )
    write_csv(work / "witness.csv", rows)
    write_csv(work / "amide.csv", amides)


if __name__ == "__main__":
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("work", type=Path)
    a.add_argument("parent", type=Path)
    a.add_argument("fixed", type=Path)
    a.add_argument("parser", type=Path)
    opts = a.parse_args()
    probe(
        opts.work.resolve(),
        opts.parent.resolve(),
        opts.fixed.resolve(),
        opts.parser.resolve(),
    )
