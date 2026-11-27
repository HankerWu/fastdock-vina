"""Re-evaluate every saved pose with an independent bound SDF and Meeko SMILES indices."""

import argparse
import csv
import hashlib
import json
import re
import sys
from importlib.metadata import version
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from scipy.spatial.transform import Rotation

from prepare_inputs import BENCH, SYSTEMS

sys.path.insert(0, str(BENCH / "pair_gradient/pose"))
from paired import write
from sources import crystal_contacts

INSTANCES = {"3P6D": 133, "3P6E": 133, "3H82": 351, "3F1O": 1}


def records(path):
    with path.open() as fh:
        return list(csv.DictReader(fh))


def xyz(line):
    return np.array([float(line[k : k + 8]) for k in (30, 38, 46)])


def independent_reference(work, code, systems=SYSTEMS, instances=INSTANCES):
    import urllib.request

    _target, component, chain = systems[code]
    instance = instances[code]
    url = (
        f"https://models.rcsb.org/v1/{code.lower()}/ligand?auth_asym_id={chain}"
        f"&auth_seq_id={instance}&encoding=sdf&copy_all_categories=false"
    )
    raw = work / "raw"
    sdf = raw / f"{code}-bound.sdf"
    if not sdf.exists():
        sdf.write_bytes(urllib.request.urlopen(url, timeout=60).read())
    bound = Chem.MolFromMolFile(str(sdf))
    Chem.AssignStereochemistryFrom3D(bound)
    lines = (work / code / "ligand.pdbqt").read_text().splitlines()
    smiles = next(
        line.split()[2] for line in lines if line.startswith("REMARK SMILES ")
    )
    mol = Chem.MolFromSmiles(smiles)
    serial = {}
    for line in lines:
        if line.startswith("REMARK SMILES IDX "):
            numbers = list(map(int, line.split()[3:]))
            serial.update(zip(numbers[::2], numbers[1::2]))
    matches = bound.GetSubstructMatches(
        mol, uniquify=False, useChirality=True, maxMatches=10001
    )
    assert (
        0 < len(matches) < 10001
        and len(serial) == mol.GetNumAtoms() == bound.GetNumAtoms()
    )
    pdb = raw / f"{code}.pdb"
    original = pdb.read_text().splitlines()
    lig = [
        line
        for line in original
        if line.startswith("HETATM")
        and line[21] == chain
        and line[17:20] == component
        and int(line[22:26]) == instance
        and line[76:78].strip() not in {"H", "D"}
    ]
    ca = {
        line[22:27].strip() + ":" + line[17:20]: xyz(line)
        for line in original
        if line.startswith("ATOM  ")
        and line[21] == chain
        and line[16] in " A"
        and line[12:16].strip() == "CA"
    }
    # ModelServer returns one conformation. Transfer its *names* via the deposited
    # atom site, then reconstruct each full alternative with its original occupancy.
    names = []
    for atom, pos in zip(bound.GetAtoms(), bound.GetConformer().GetPositions()):
        ids = {
            line[12:16].strip()
            for line in lig
            if line[76:78].strip().upper() == atom.GetSymbol().upper()
            and np.linalg.norm(xyz(line) - pos) < 0.002
        }
        if len(ids) != 1:
            raise ValueError(
                f"Ambiguous ModelServer/PDB coordinate correspondence {code}: {ids}"
            )
        names.append(ids.pop())
    assert len(set(names)) == len(names)
    ref, labels, occ = [], [], []
    alts = sorted({line[16] for line in lig} - {" "}) or [" "]
    for alt in alts:
        coords = {
            line[12:16].strip(): xyz(line) for line in lig if line[16] in (" ", alt)
        }
        occupancy = {float(line[54:60]) for line in lig if line[16] == alt}
        assert len(occupancy) == 1
        occ.append(occupancy.pop())
        ref.extend(np.asarray([coords[n] for n in names])[np.asarray(matches)])
        labels.extend([alt.strip() or "."] * len(matches))
    for row in original:
        if row.startswith("REMARK   2 RESOLUTION."):
            resolution = float(row.split()[3])
            break
    info = dict(
        source=url,
        sha256=hashlib.sha256(sdf.read_bytes()).hexdigest(),
        altloc=alts,
        occupancies=occ,
        resolution_A=resolution,
        heavy_atoms=len(names),
        stereochemical_matches=len(matches),
        chain=chain,
        formal_charge=Chem.GetFormalCharge(mol),
        **crystal_contacts(gemmi.read_structure(str(pdb)), chain, component),
    )
    center = np.array(ref[0]).mean(0)
    return (
        np.asarray(ref),
        ca,
        center,
        [serial[i] for i in range(1, len(serial) + 1)],
        labels,
        info,
    )


def main(work, results, instances=INSTANCES, systems=SYSTEMS):
    sources = {
        code: independent_reference(work, code, systems, instances)
        for code in instances
    }
    per = records(results / "cells.csv")
    modes = {
        (r["key"], int(r["mode"])): float(r["rmsd_A"])
        for r in records(results / "modes.csv")
    }
    checked, max_difference = [], 0
    for cell in per:
        refs, ca_q, cent_q, serial, labels, _ = sources[cell["query"]]
        _, ca_r, cent_r, *_ = sources[cell["receiver"]]
        common = sorted(
            k
            for k in ca_q.keys() & ca_r.keys()
            if np.linalg.norm(ca_q[k] - cent_q) < 14
            and np.linalg.norm(ca_r[k] - cent_r) < 14
        )
        x, y = [np.array([ca[k] for k in common]) for ca in (ca_q, ca_r)]
        rot, _ = Rotation.align_vectors(y - y.mean(0), x - x.mean(0))
        aligned = rot.apply((refs - x.mean(0)).reshape(-1, 3)).reshape(
            refs.shape
        ) + y.mean(0)
        key = cell["key"]
        file = results / "runs" / (key + ".pdbqt")
        assert hashlib.sha256(file.read_bytes()).hexdigest() == cell["output_sha256"]
        for i, block in enumerate(file.read_text().split("ENDMDL")[:-1], 1):
            coords = {
                int(l[6:11]): xyz(l)
                for l in block.splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
            }
            pred = np.array([coords[j] for j in serial])
            values = np.linalg.norm(pred - aligned, axis=(1, 2)) / np.sqrt(len(serial))
            rms = values.min()
            err = abs(rms - modes[key, i])
            if err > 1e-6:
                raise ValueError(f"{key}:{i} coordinate mismatch {err} A")
            max_difference = max(max_difference, err)
            checked.append(
                dict(
                    key=key,
                    mode=i,
                    independent_A=rms,
                    best_alt=labels[np.argmin(values)],
                    alt_A=min(v for v, l in zip(values, labels) if l in (".", "A")),
                    alt_B=min(
                        (v for v, l in zip(values, labels) if l == "B"), default=""
                    ),
                )
            )
        log = (results / "runs" / (key + ".counts")).read_text()
        for name in (
            "hunt_calls",
            "hunt_evals",
            "second_calls",
            "second_evals",
            "explicit_inputs",
            "explicit_calls",
            "explicit_evals",
        ):
            ns = re.findall(r"\b" + name + r"=(\d+)\b", log)
            assert len(ns) == (1 if name.startswith("explicit") else 4) and sum(
                map(int, ns)
            ) == int(cell[name])
    software = {
        v: version(v)
        for v in ("rdkit", "meeko", "gemmi", "numpy", "scipy", "openmm", "pdbfixer")
    }
    (results / "verification.json").write_text(
        json.dumps(
            dict(
                cells=len(per),
                independent_modes=len(checked),
                max_difference_A=max_difference,
                six_work_counters_and_input_counts_verified=True,
                inputs={k: v[-1] for k, v in sources.items()},
                software=software,
            ),
            indent=2,
        )
        + "\n"
    )
    write(results / "independent_modes.csv", checked)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("results", type=Path)
    a = parser.parse_args()
    main(a.work, a.results)
