"""Reparse coordinates via the independent receptor fit and PDBQT canonical SMILES order."""

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from urllib.request import urlopen

from panel import HERE, PANELS
from rdkit import Chem

sys.path.insert(0, str(HERE.parent / "paired_crossdocking"))
import independent_cross as check
import numpy as np


def instance(work, code, comp):
    lines = (work / "raw" / (code + ".pdb")).read_text().splitlines()
    coords = {
        (l[17:20].strip(), l[22:26], l[26]): check.xyz(l)
        for l in lines
        if l.startswith("ATOM")
        and l[21] == "A"
        and l[16] in " A"
        and l[12:16].strip() == "CA"
    }
    observed = [
        l
        for l in lines
        if l.startswith("HETATM")
        and l[21] == "A"
        and l[17:20] == comp
        and l[16] in " A"
        and l[76:78].strip() not in ("H", "D")
    ]
    positions = {int(l[22:26]) for l in observed}
    if len(positions) != 1:
        raise ValueError(f"Ambiguous instance {code}")
    uri = f"https://models.rcsb.org/v1/{code.lower()}/ligand?auth_asym_id=A&auth_seq_id={positions.pop()}&encoding=sdf"
    sdf = work / "raw" / (code + "_instance.sdf")
    if not sdf.exists():
        sdf.write_bytes(urlopen(uri, timeout=60).read())
    bound = Chem.MolFromMolFile(str(sdf), removeHs=True)
    if bound is None:
        raise ValueError(f"Invalid instance SDF for {code}")
    Chem.AssignStereochemistryFrom3D(bound)
    query_lines = (work / code / "ligand.pdbqt").read_text().splitlines()
    text = next(
        l.split(maxsplit=2)[2]
        for l in query_lines
        if l.startswith("REMARK SMILES ") and "IDX" not in l
    )
    query = Chem.MolFromSmiles(text)
    idx = {}
    for l in query_lines:
        if l.startswith("REMARK SMILES IDX"):
            numbers = list(map(int, l.split()[3:]))
            idx.update(zip(numbers[::2], numbers[1::2]))
    prepared = {
        int(l[6:11]): check.xyz(l)
        for l in query_lines
        if l.startswith(("ATOM  ", "HETATM"))
    }
    conf = Chem.Conformer(query.GetNumAtoms())
    for i in range(query.GetNumAtoms()):
        conf.SetAtomPosition(i, prepared[idx[i + 1]])
    query.AddConformer(conf)
    # RDKit/Meeko may serialize a symmetry-related ring traversal. Resolve
    # declared stereocentres in that actual serial order using the input geometry.
    declared = {
        a.GetIdx()
        for a in query.GetAtoms()
        if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    }
    Chem.AssignStereochemistryFrom3D(query)
    for atom in query.GetAtoms():
        if atom.GetIdx() not in declared:
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    match = bound.GetSubstructMatches(
        query, useChirality=True, uniquify=False, maxMatches=10000
    )
    if not match or len(match) == 10000 or query.GetNumAtoms() != bound.GetNumAtoms():
        raise ValueError(f"Unmatched chemical graph {code}")
    mapped = [bound.GetConformer().GetPositions()[list(m)] for m in match]
    return (
        coords,
        np.mean([check.xyz(l) for l in observed], axis=0),
        mapped,
        [idx[i + 1] for i in range(query.GetNumAtoms())],
        mapped,
    ), (uri, hashlib.sha256(sdf.read_bytes()).hexdigest())


def verify(work, result):
    sources, data = {}, {}
    for group in PANELS.values():
        for code, comp in group.items():
            data[code], (uri, digest) = instance(work, code, comp)
            sources[uri] = digest
    (work / "instance_sources.json").write_text(json.dumps(sources, indent=2) + "\n")
    max_rms = max_align = 0
    count = 0
    key = lambda row: (
        row["query"],
        row["receiver"],
        row["cap"],
        row["seed"],
        row["engine"],
        int(row["model"]),
    )
    with (result / "poses.csv").open() as stream:
        expected = {key(row): row for row in csv.DictReader(stream)}
    with (work / "geometry.csv").open() as stream:
        geom = {key(row): row for row in csv.DictReader(stream)}
    for path in sorted((work / "runs").glob("*.pdbqt")):
        src, rec, cap, seed, engine = path.stem.split("_")
        refs, _, n, ca_rms = check.fitted(data[src], data[rec])
        serials = data[src][3]
        for block in path.read_text().split("ENDMDL"):
            lines = block.strip().splitlines()
            if not lines:
                continue
            rank = int(lines[0].split()[1])
            item = (src, rec, cap, seed, engine, rank)
            published = expected.pop(item)
            xyz = {
                int(l[6:11]): check.xyz(l)
                for l in lines
                if l.startswith(("ATOM", "HETATM"))
            }
            measured = min(
                np.linalg.norm(ref - np.array([xyz[s] for s in serials]))
                / math.sqrt(len(serials))
                for ref in refs
            )
            max_rms = max(max_rms, abs(float(published["rmsd"]) - measured))
            if abs(float(published["rmsd"]) - measured) > 1e-7:
                raise ValueError(
                    f"{item}: independent RMSD {measured}, reported {published['rmsd']}"
                )
            max_align = max(max_align, abs(float(geom[item]["ca_rmsd"]) - ca_rms))
            assert int(geom[item]["n_ca"]) == n
            count += 1
    assert not expected and max_rms < 1e-7 and max_align < 1e-8
    print(
        f"Verified {count} poses. Maximum RMSD discrepancy {max_rms:.3g} Å, maximum receptor fit discrepancy {max_align:.3g} Å."
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("result", type=Path)
    args = p.parse_args()
    verify(args.work, args.result)
