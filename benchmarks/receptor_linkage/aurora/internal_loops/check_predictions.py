"""Independent checks of stereochemistry, fixed scores and archived blind poses."""

import argparse
import csv
import hashlib
import json
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import gemmi
import numpy as np
from rdkit.Chem import rdMolAlign
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import frame_check as frame

NEW = set(range(281, 292)) | {303, 304, 305}
LINKS = sorted(NEW | {280, 302})
DRUGS = {"danusertib": ("2J50", "627"), "at9283": ("2W1G", "L0G")}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atoms(path):
    xyz = {}
    for line in path.read_text().splitlines():
        if (
            line.startswith("ATOM  ")
            and line[21] == "A"
            and line[16] in " A"
            and line[76:78].strip() not in ("H", "D")
        ):
            key = int(line[22:26]), line[17:20], line[12:16].strip()
            if key in xyz:
                raise ValueError(f"Duplicate {key}")
            xyz[key] = np.array([float(line[k : k + 8]) for k in (30, 38, 46)])
    return xyz


def centres(coords):
    by_res = defaultdict(dict)
    for (i, res, atom), pos in coords.items():
        by_res[i, res][atom] = pos
    signed = {}
    for (i, res), xyz in by_res.items():
        for site, names in [("CA", ("N", "C", "CB")), ("CB", ("CA", "OG1", "CG2"))]:
            if site in xyz and all(a in xyz for a in names):
                signed[i, res, site] = float(
                    np.linalg.det([xyz[a] - xyz[site] for a in names])
                )
    return signed


def geometry(path, raw):
    from openmm.app import PDBFile

    xyz = atoms(path)
    hand = centres(xyz)
    observed = centres(raw)
    ref = {
        site: int(np.sign(np.median([v for k, v in observed.items() if k[-1] == site])))
        for site in ("CA", "CB")
    }
    bad = [
        f"{i}-{name}-{site}"
        for (i, name, site), vol in hand.items()
        if np.sign(vol) != ref[site]
    ]
    by_index = {(i, a): v for (i, _, a), v in xyz.items()}
    omega = {}
    for i in LINKS:
        p = [
            by_index[j, a]
            for j, a in [(i, "CA"), (i, "C"), (i + 1, "N"), (i + 1, "CA")]
        ]
        axis = p[2] - p[1]
        axis /= np.linalg.norm(axis)
        b, c = p[0] - p[1], p[-1] - p[2]
        x, y = b - np.dot(b, axis) * axis, c - np.dot(c, axis) * axis
        omega[f"{i}-{i + 1}"] = float(
            np.rad2deg(np.arctan2(np.dot(np.cross(axis, x), y), np.dot(x, y)))
        )
    # Standard residue bonds define graph exclusions independently of length.
    pdb = PDBFile(str(path))
    network = defaultdict(set)
    for a, b in pdb.topology.bonds():
        i = (int(a.residue.id), a.name)
        j = (int(b.residue.id), b.name)
        if i in by_index and j in by_index:
            network[i].add(j)
            network[j].add(i)
    keys = list(by_index)
    tree = cKDTree([by_index[k] for k in keys])
    below = []
    for h, l in tree.query_pairs(2.0):
        a, b = keys[h], keys[l]
        if a[0] not in NEW and b[0] not in NEW:
            continue
        excluded = {a}
        for _ in range(3):
            excluded |= {v for u in excluded for v in network[u]}
        if b not in excluded:
            below.append((float(np.linalg.norm(by_index[a] - by_index[b])), a, b))
    gap_lengths = [
        np.linalg.norm(by_index[i, "C"] - by_index[i + 1, "N"]) for i in LINKS
    ]
    return {
        "sha256": digest(path),
        "centres_checked_new": sum(i in NEW for i, _, site in hand),
        "native_sign_N_C_CB_at_CA": ref["CA"],
        "native_sign_CA_OG1_CG2_at_CB": ref["CB"],
        "wrong_centres_in_chain": bad,
        "new_centre_determinants_A3": {
            f"{i}-{n}-{a}": v for (i, n, a), v in hand.items() if i in NEW
        },
        "leu_CG": "Two terminal methyl groups CD1 and CD2 are equivalent; no additional stereocentre.",
        "peptide_omega_deg": omega,
        "nonpro_cis_like_abs_omega_lt_30": [
            k for k, v in omega.items() if abs(v) < 30 and k != "281-282"
        ],
        "max_abs_deviation_from_cis_or_trans_deg": max(
            min(abs(v), 180 - abs(v)) for v in omega.values()
        ),
        "peptide_CN_A": {
            f"{i}-{i + 1}": float(v) for i, v in zip(LINKS, gap_lengths, strict=True)
        },
        "peptide_CN_range_A": [float(min(gap_lengths)), float(max(gap_lengths))],
        "heavy_pairs_lt_2A_graph_separation_gt_3_involving_new_residues": sorted(below),
    }


def template_check(work, folder, raw):
    from openmm.app import PDBFile, forcefield

    pdb = folder / "repaired.pdb"
    report = geometry(pdb, raw)
    modeled = atoms(pdb)
    source = atoms(work / "prepared/2X6E-protein.pdb")
    old = {k: v for k, v in modeled.items() if k[0] not in NEW}
    assert set(old) == set(source)
    max_move = max(np.linalg.norm(v - source[k]) for k, v in old.items())
    report["historical_max_heavy_displacement_A"] = float(max_move)
    assert max_move < 0.0011
    xml = ET.parse(
        Path(forcefield.__file__).parent / "data/amber14/protein.ff14SB.xml"
    ).getroot()
    types = {
        (r.get("name"), a.get("name")): a.get("type")
        for r in xml.find("Residues")
        for a in r.findall("Atom")
    }
    r0 = {
        tuple(sorted([b.get("type1"), b.get("type2")])): 10 * float(b.get("length"))
        for b in xml.find("HarmonicBondForce")
    }
    deviations = []
    for a, b in PDBFile(str(pdb)).topology.bonds():

        def ident(atom):
            return (int(atom.residue.id), atom.residue.name, atom.name)

        ka, kb = ident(a), ident(b)
        if not (ka[0] in NEW or kb[0] in NEW) or ka not in modeled or kb not in modeled:
            continue
        ts = tuple(
            sorted(
                types.get(
                    (k[1], k[2]), types.get(("HIE", k[2])) if k[1] == "HIS" else None
                )
                for k in (ka, kb)
            )
        )
        deviations.append(
            abs(float(np.linalg.norm(modeled[ka] - modeled[kb])) - r0[ts])
        )
    report.update(
        all_new_heavy_bonds_checked=len(deviations),
        maximum_ff14sb_bond_deviation_A=max(deviations),
    )
    if not (folder / "receptor.pdbqt").exists():
        report["typed_for_docking"] = False
        return report
    historical = (work / "prepared/2X6E.pdbqt").read_bytes()
    complete = (folder / "receptor.pdbqt").read_bytes()
    assert complete.startswith(historical)
    added = complete[len(historical) :].decode().splitlines()
    typed = {
        l[12:]
        for l in (folder / "typed.pdbqt").read_text().splitlines()
        if l.startswith("ATOM") and int(l[22:26]) in NEW
    }
    assert len(added) == len(typed) and typed == {l[12:] for l in added}
    observed = {k for k in modeled if k[0] in NEW}
    actual = {
        (int(l[22:26]), l[17:20], l[12:16].strip())
        for l in added
        if l.split()[-1] not in ["H", "HD"]
    }
    assert observed == actual
    assert not report["wrong_centres_in_chain"] and max(deviations) <= 0.15
    assert not report["heavy_pairs_lt_2A_graph_separation_gt_3_involving_new_residues"]
    assert all(1.2 <= v <= 1.47 for v in report["peptide_CN_A"].values())
    assert all(abs(v) >= 145 for v in report["peptide_omega_deg"].values())
    report.update(
        typed_for_docking=True,
        unchanged_prior_record_count=len(historical.splitlines()),
        new_typed_records=len(added),
        new_heavy_records=len(actual),
        receptor_sha256=digest(folder / "receptor.pdbqt"),
    )
    return report


def refmol(work, drug):
    code, component = DRUGS[drug]
    original = gemmi.read_structure(str(work / "raw/2X6E.pdb"))
    original.merge_chain_parts()
    receiver = original[0]["A"]
    path = work / (code + ".pdb")
    crystal = gemmi.read_structure(str(path))[0]["A"]
    fixed, moving = frame.ca(receiver), frame.ca(crystal)
    box = [
        a.pos
        for r in receiver
        if r.name == "YM4"
        for a in r
        if not a.element.is_hydrogen and a.altloc in ("\x00", "A")
    ]
    keys = [
        k
        for k in sorted(fixed.keys() & moving.keys())
        if min(fixed[k].dist(p) for p in box) <= 12
    ]
    fit = gemmi.superpose_positions([fixed[k] for k in keys], [moving[k] for k in keys])
    reference = frame.reference_graph(path, component, "A")
    conf = reference.GetConformer()
    for i, pos in enumerate(conf.GetPositions()):
        conf.SetAtomPosition(i, tuple(fit.transform.apply(gemmi.Position(*pos))))
    return reference, len(keys)


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def audit_blind(work, result_dirs):
    output = {}
    for directory in result_dirs:
        if not (directory / "cells.csv").exists():
            continue
        cells = rows(directory / "cells.csv")
        published = {
            (r["arm"], r["drug"], r["seed"], int(r["rank"])): r
            for r in rows(directory / "modes.csv")
        }
        mode_errors, cell_errors, binding_errors = [], [], []
        for drug in DRUGS:
            ref, _ = refmol(work, drug)
            mol, serials = frame.read_query(work, drug)
            for row in (r for r in cells if r["drug"] == drug):
                arm, seed = row["arm"], row["seed"]
                path = directory / "blind" / f"{arm}_{drug}_{seed}.pdbqt"
                if digest(path) != row["pool_sha256"]:
                    raise ValueError(f"Pose archive mismatch {path}")
                coords = frame.positions(path)
                scores = [
                    float(l.split()[3])
                    for l in path.read_text().splitlines()
                    if l.startswith("REMARK VINA RESULT:")
                ]
                err = []
                for rank, (model, score) in enumerate(
                    zip(coords, scores, strict=True), 1
                ):
                    for i, s in enumerate(serials):
                        mol.GetConformer().SetAtomPosition(i, model[s])
                    value = rdMolAlign.CalcRMS(
                        mol, ref, symmetrizeConjugatedTerminalGroups=False
                    )
                    err.append(value)
                    pub = published[arm, drug, seed, rank]
                    mode_errors.append(abs(value - float(pub["rmsd_A"])))
                    binding_errors.append(abs(score - float(pub["score"])))
                assert len(err) == int(row["n"])
                derived = (err[0], min(err[:5]), min(err), err[0] - min(err))
                cell_errors.extend(
                    abs(v - float(row[k]))
                    for v, k in zip(
                        derived, ["top1", "best5", "oracle", "regret"], strict=True
                    )
                )
        assert (
            len(published) == len(mode_errors)
            and max(mode_errors + binding_errors + cell_errors) < 1e-7
        )
        output[directory.name] = {
            "cells": len(cells),
            "modes": len(mode_errors),
            "max_mode_rmsd_difference_A": max(mode_errors),
            "max_score_difference": max(binding_errors),
            "max_cell_metric_difference_A": max(cell_errors),
        }
    return output


def fixed(work, directory):
    import re

    source = directory / "fixed.csv"
    if not source.exists():
        return None
    records = rows(source)
    gaps = {}
    for drug in DRUGS:
        ref, _ = refmol(work, drug)
        mol, serials = frame.read_query(work, drug)
        for arm in sorted({r["arm"] for r in records}):
            subset = [r for r in records if r["drug"] == drug and r["arm"] == arm]
            scores = {}
            for number, r in enumerate(subset, 1):
                filename = (
                    work / "coordinates" / f"2X6E_{drug}_2X6E_pocket_{number}.pdbqt"
                )
                assert digest(filename) == r["input_sha256"]
                # Coordinate witnesses have no MODEL/ENDMDL wrappers.
                coordinates = {
                    int(s[6:11]): [float(s[k : k + 8]) for k in (30, 38, 46)]
                    for s in filename.read_text().splitlines()
                    if s.startswith(("ATOM", "HETATM"))
                }
                for i, s in enumerate(serials):
                    mol.GetConformer().SetAtomPosition(i, coordinates[s])
                rmsd = rdMolAlign.CalcRMS(
                    mol, ref, symmetrizeConjugatedTerminalGroups=False
                )
                assert abs(rmsd - float(r["rmsd_A"])) < 1e-7
                log = (
                    directory / "fixed" / f"{drug}_{arm}_{r['pose']}.log"
                ).read_text()
                score = float(
                    re.search(
                        r"Estimated Free Energy of Binding\s*:\s*([-+.\d]+)", log
                    )[1]
                )
                assert score == float(r["score"])
                scores[r["pose"]] = score
            assert len(scores) == 5
            gaps[drug, arm] = {
                p: round(s - scores["witness"], 3)
                for p, s in scores.items()
                if p != "witness"
            }
            for r in subset:
                assert (
                    abs(
                        scores[r["pose"]]
                        - scores["witness"]
                        - float(r["decoy_minus_witness"])
                    )
                    < 1e-7
                )
    return {
        "decoy_minus_witness": {f"{d}_{a}": v for (d, a), v in gaps.items()},
        "change_in_gap_relative_to_dry": {
            f"{drug}_{arm}": {
                p: round(value - gaps[drug, "dry"][p], 3) for p, value in values.items()
            }
            for (drug, arm), values in gaps.items()
            if arm != "dry" and (drug, "dry") in gaps
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument(
        "--results",
        nargs="*",
        type=Path,
        help="directories containing cells/modes and/or fixed.csv",
    )
    args = parser.parse_args()
    work = args.work
    dirs = args.results or [work / "base"]
    raw = atoms(work / "raw/2X6E.pdb")
    models = {}
    for folder in sorted(work.glob("loop*")):
        if folder.is_dir():
            for path in sorted(folder.glob("*.pdb")):
                models[str(path.relative_to(work))] = geometry(path, raw)
    templates = {
        str(folder.relative_to(work)): template_check(work, folder, raw)
        for folder in sorted(work.glob("template*/trans_model/*"))
        if (folder / "repaired.pdb").exists()
    }
    report = {
        "template_controls": templates,
        "models": models,
        "blind": audit_blind(work, dirs),
        "fixed": {p.name: fixed(work, p) for p in dirs},
    }
    print(json.dumps(report, indent=2, allow_nan=False))
