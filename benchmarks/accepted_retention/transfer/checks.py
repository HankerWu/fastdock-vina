"""Reconstruct structures from raw PDB/CCD, and reconcile counts and shared modes."""

import argparse
import csv
import hashlib
import json
from collections import Counter
from functools import cache
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation
from inputs import CASES, reference as engine
from audit import coordinates, reference


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def analyse(work):
    info = {p: reference(work, p, ccd, ch) for p, (_, ccd, ch) in CASES.items()}

    @cache
    def frame(q, r):
        qa, qc, mol = info[q]
        ra, rc, _ = info[r]
        keys = sorted(
            k
            for k in qa.keys() & ra.keys()
            if np.linalg.norm(qa[k] - qc) < 14 and np.linalg.norm(ra[k] - rc) < 14
        )
        x = np.stack([qa[k] for k in keys])
        y = np.stack([ra[k] for k in keys])
        rot, _ = Rotation.align_vectors(y - y.mean(0), x - x.mean(0))
        return mol, rot.apply(mol.GetConformer().GetPositions() - x.mean(0)) + y.mean(0)

    atom_refs = {}
    for q in CASES:
        mol, xyz, tags = coordinates(
            (work / q / "ligand.pdbqt").read_text().splitlines()
        )
        bonds = np.array(
            [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
        )
        atom_refs[q] = (
            tags,
            bonds,
            np.linalg.norm(xyz[bonds[:, 0]] - xyz[bonds[:, 1]], axis=1),
        )
    deviations = []
    bond_error = 0
    inverted = 0
    fingerprints = {}
    modes = {}
    for file in ("modes.csv", "release_modes.csv", "time_modes.csv"):
        if (work / file).exists():
            for row in read(work / file):
                modes[row["key"], int(row["mode"])] = (
                    float(row["rmsd_A"]),
                    float(row["score"]),
                )
    for key in sorted({k for k, _ in modes}):
        q, r, *_ = key.split("_")
        mol, ref = frame(q, r)
        fingerprints[key] = []
        initial_tags, bonds, lengths = atom_refs[q]
        for n, block in enumerate(
            (work / "runs" / (key + ".pdbqt")).read_text().split("ENDMDL")[:-1], 1
        ):
            ls = block.splitlines()
            probe, xyz, tags = coordinates(ls)
            maps = mol.GetSubstructMatches(
                probe, uniquify=False, useChirality=True, maxMatches=10001
            )
            assert 0 < len(maps) < 10001
            dist = min(
                np.linalg.norm(xyz - ref[np.array(m)]) / np.sqrt(len(m)) for m in maps
            )
            refdist, score = modes.pop((key, n))
            assert (
                float(
                    next(l for l in ls if l.startswith("REMARK VINA RESULT")).split()[3]
                )
                == score
            )
            deviations.append(abs(dist - refdist))
            bond_error = max(
                bond_error,
                float(
                    np.max(
                        abs(
                            lengths
                            - np.linalg.norm(
                                xyz[bonds[:, 0]] - xyz[bonds[:, 1]], axis=1
                            )
                        )
                    )
                ),
            )
            inverted += int(
                any(
                    t != "CHI_UNSPECIFIED" and tags[i] != t
                    for i, t in initial_tags.items()
                )
            )
            fingerprints[key].append(hashlib.sha256(xyz.tobytes()).hexdigest())
    assert not modes and max(deviations) < 1e-6
    original = read(work / "cells.csv")
    groups = {}
    for row in original:
        arm = row["arm"]
        key = row["key"]
        log = work / "runs" / (key + (".counts" if arm == "upstream" else ".log"))
        fromlog = engine.counts(log.read_text())
        assert all(int(row[k]) == int(v) for k, v in fromlog.items())
        groups.setdefault(key.rsplit("_", 1)[0], {})[arm] = row
    pairs = []
    for _, by in sorted(groups.items()):
        b = by["parent"]
        baseline = fingerprints[b["key"]]
        for arm in ("gated", "unconditional", "upstream"):
            a = by[arm]
            got = fingerprints[a["key"]]
            same_search = all(
                a[k] == b[k]
                for k in ("hunt_calls", "hunt_evals", "second_calls", "second_evals")
            )
            if arm != "upstream":
                assert same_search and a["explicit_inputs"] == b["explicit_inputs"]
                assert int(a["replaced"]) or baseline == got
            pairs.append(
                dict(
                    key=a["key"],
                    search_counts_equal=same_search,
                    shared_modes=sum((Counter(baseline) & Counter(got)).values()),
                    parent_modes=len(baseline),
                    arm_modes=len(got),
                    identical_file=(work / "runs" / (a["key"] + ".pdbqt")).read_bytes()
                    == (work / "runs" / (b["key"] + ".pdbqt")).read_bytes(),
                    top1_coordinates_equal=baseline[0] == got[0],
                    top1_delta_A=float(a["top1_A"]) - float(b["top1_A"]),
                    oracle_delta_A=float(a["oracle_A"]) - float(b["oracle_A"]),
                )
            )
    engine.write(work / "paired_changes.csv", pairs)
    return dict(
        poses_verified=len(deviations),
        max_rmsd_difference_A=max(deviations),
        max_covalent_bond_length_change_A=bond_error,
        poses_with_stereocenter_tag_change=inverted,
        source_count_log_checks=len(original),
        paired_comparisons=len(pairs),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    a = p.parse_args()
    (a.work / "audit.json").write_text(json.dumps(analyse(a.work), indent=2) + "\n")
