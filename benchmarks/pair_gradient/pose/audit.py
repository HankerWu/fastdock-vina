"""Check every pose independently from deposited chemistry and backbone coordinates."""

import argparse
import csv
import json
import sys
from functools import cache
from pathlib import Path

import gemmi
import numpy as np
from Bio import Align
from prepare import CASES, ROOT
from rdkit import Chem
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(ROOT / "benchmarks/cubic_grid/human_panel"))
from verify import reference


def coordinates(lines):
    smiles = next(
        s.split(maxsplit=2)[2]
        for s in lines
        if s.startswith("REMARK SMILES ") and not s.startswith("REMARK SMILES IDX")
    )
    mol = Chem.MolFromSmiles(smiles)
    mapping = {}
    for s in lines:
        if s.startswith("REMARK SMILES IDX"):
            nums = [int(x) for x in s.split()[3:]]
            mapping.update(zip(nums[::2], nums[1::2], strict=True))
    positions = {
        int(s[6:11]): [float(s[j : j + 8]) for j in (30, 38, 46)]
        for s in lines
        if s.startswith(("ATOM  ", "HETATM"))
    }
    xyz = np.array([positions[mapping[j + 1]] for j in range(mol.GetNumAtoms())])
    conf = Chem.Conformer(len(xyz))
    for i, x in enumerate(xyz):
        conf.SetAtomPosition(i, x)
    mol.AddConformer(conf)
    stereo = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(stereo, replaceExistingTags=True)
    tags = {a.GetIdx(): str(a.GetChiralTag()) for a in stereo.GetAtoms()}
    return mol, xyz, tags


def align_renin_numbering(work, info, cases):
    canonical = gemmi.read_structure(str(work / "raw/1RNE.pdb"))
    canonical.setup_entities()
    seq = gemmi.one_letter_code(canonical.entities[0].full_sequence)
    aligner = Align.PairwiseAligner(
        mode="global",
        match_score=2,
        mismatch_score=-2,
        open_gap_score=-8,
        extend_gap_score=-0.5,
    )
    for pdb, (target, _, _) in cases.items():
        if target != "REN":
            continue
        ca, center, mol = info[pdb]
        keys = list(ca)
        observed = gemmi.one_letter_code([k[-1] for k in keys])
        aligned = aligner.align(seq, observed)[0]
        positions = {}
        for (a, b), (x, y) in zip(*aligned.aligned):
            for i, j in zip(range(a, b), range(x, y), strict=True):
                if seq[i] == observed[j]:
                    positions[(i, keys[j][-1])] = ca[keys[j]]
        info[pdb] = positions, center, mol


def main(work, cases=CASES):
    info = {p: reference(work, p, c, ch) for p, (_, c, ch) in cases.items()}
    if any(target == "REN" for target, _, _ in cases.values()):
        align_renin_numbering(work, info, cases)

    @cache
    def receiver_frame(query, receptor):
        qa, qc, mol = info[query]
        ra, rc, _ = info[receptor]
        keys = sorted(
            k
            for k in qa.keys() & ra.keys()
            if np.linalg.norm(qa[k] - qc) < 14 and np.linalg.norm(ra[k] - rc) < 14
        )
        x, y = np.stack([qa[k] for k in keys]), np.stack([ra[k] for k in keys])
        rotation, _ = Rotation.align_vectors(y - y.mean(0), x - x.mean(0))
        return mol, rotation.apply(
            mol.GetConformer().GetPositions() - x.mean(0)
        ) + y.mean(0)

    input_geometry = {}
    for pdb in cases:
        mol, xyz, tags = coordinates(
            (work / pdb / "ligand.pdbqt").read_text().splitlines()
        )
        pairs = np.array(
            [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
        )
        length = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1)
        input_geometry[pdb] = pairs, length, tags

    results = {}
    for name in (
        "modes.csv",
        "matched_modes.csv",
        "shift_modes.csv",
        "release_modes.csv",
    ):
        path = work / name
        if path.exists():
            with path.open() as f:
                for r in csv.DictReader(f):
                    results[r["key"], int(r["mode"])] = float(r["rmsd_A"])
    largest = bond_error = 0.0
    checked = inversions = 0
    for key in sorted({k for k, _ in results}):
        query, receiver, *_ = key.split("_")
        mol, ref = receiver_frame(query, receiver)
        pairs, lengths, input_tags = input_geometry[query]
        for num, block in enumerate(
            (work / "runs" / (key + ".pdbqt")).read_text().split("ENDMDL")[:-1], 1
        ):
            m, xyz, tags = coordinates(block.splitlines())
            inversions += int(
                any(
                    tag != "CHI_UNSPECIFIED" and tags[i] != tag
                    for i, tag in input_tags.items()
                )
            )
            trial_length = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1)
            bond_error = max(bond_error, float(np.max(abs(trial_length - lengths))))
            maps = mol.GetSubstructMatches(
                m, uniquify=False, useChirality=True, maxMatches=10001
            )
            if not 0 < len(maps) < 10001:
                raise ValueError((key, num, "invalid stereo-sensitive correspondence"))
            dist = min(
                np.linalg.norm(ref[list(match)] - xyz) / np.sqrt(len(match))
                for match in maps
            )
            largest = max(largest, abs(dist - results.pop((key, num))))
            checked += 1
    if results or largest > 1e-5:
        raise ValueError((largest, len(results)))
    return {
        "unique_modes": checked,
        "max_rmsd_error_A": largest,
        "modes_with_input_chiral_tag_changes": inversions,
        "max_bond_length_delta_from_prepared_A": bond_error,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    args = p.parse_args()
    result = main(args.work)
    (args.work / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(result)
