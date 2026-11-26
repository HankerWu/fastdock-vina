"""Independently verify Vina radial terms, HD ownership and bounded angular penalties."""
import argparse
import csv
import json
import re
import tarfile
import tempfile
from collections import defaultdict
from functools import cache
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from evaluate import INPUTS, UPSTREAM


def atoms(text):
    return {int(l[6:11]): (l[77:].strip(), l[17:27], np.array([float(l[k:k+8]) for k in (30,38,46)]))
            for l in text.splitlines() if l.startswith(("ATOM  ", "HETATM"))}


def mapping(text, label):
    result = []
    for line in text.splitlines():
        if line.startswith(label):
            values = list(map(int,line[len(label):].split()))
            result.extend(zip(values[::2],values[1::2],strict=True))
    return result


@cache
def read_inputs(folder, pdb, is_receptor):
    text = (folder / pdb / ("receptor.pdbqt" if is_receptor else "ligand.pdbqt")).read_text()
    a = atoms(text)
    h = defaultdict(list)
    if is_receptor:
        residues = defaultdict(list)
        for i,(element, residue, xyz) in a.items():
            if element not in ("H", "HD"):
                residues[residue].append(i)
        for i,(element, residue, xyz) in a.items():
            if element != "HD":
                continue
            distances = sorted((np.linalg.norm(xyz-a[j][2]),j) for j in residues[residue])
            if not distances or distances[0][0] > 1.5 or distances[1][0] - distances[0][0] < 0.1:
                raise ValueError((pdb, i, distances[:2]))
            h[distances[0][1]].append(i)
    else:
        smiles_to_serial = dict(mapping(text, "REMARK SMILES IDX "))
        for parent, serial in mapping(text, "REMARK H PARENT "):
            h[smiles_to_serial[parent]].append(serial)
    assert {i for v in h.values() for i in v} == {i for i in a if a[i][0]=="HD"}
    ids = [i for i in a if a[i][0] in ("N","NA","O","OA")]
    return a, h, ids, cKDTree([a[i][2] for i in ids])


def quality(d, h, a):
    # Direct angle at H. 180 degrees is ideal, 90 or below gives zero.
    bonds = (d-h, a-h)
    dot = np.dot(*bonds) / (np.linalg.norm(bonds[0])*np.linalg.norm(bonds[1]))
    return max(0., -dot)**2


def audit(output):
    with (output / "poses.csv").open() as f:
        table = list(csv.DictReader(f))
    groups = defaultdict(list)
    for row in table:
        groups[tuple(row[k] for k in ("protocol","design","query","receiver","seed"))].append(row)
    max_radial, max_penalty, count = 0., 0., 0
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for protocol, archive in INPUTS.items():
            with tarfile.open(archive) as t:
                t.extractall(root/protocol,filter="data")
        for protocol, design in sorted({k[:2] for k in groups}):
            source = UPSTREAM / protocol
            if design != "fixed":
                source /= design
            with tarfile.open(source/"raw.tar.xz") as t:
                for (family, batch, q, r, seed), values in groups.items():
                    if (protocol, design) != (family, batch):
                        continue
                    lig, lh, lids, _ = read_inputs(root/family,q,False)
                    rec, rh, rids, tree = read_inputs(root/family,r,True)
                    raw = t.extractfile(f"{q}-{r}-{seed}-native.pdbqt").read().decode()
                    blocks = re.findall(r"^MODEL[^\n]*\n(.*?)^ENDMDL$",raw,re.M | re.S)
                    for block, value in zip(blocks,values,strict=True):
                        xyz = atoms(block)
                        native_contact = 0
                        radial, penalty = 0., 0.
                        for i in lids:
                            li = xyz[i][2]
                            for j0 in tree.query_ball_point(li,3.6):
                                j = rids[j0]
                                rj = rec[j][2]
                                orient = []
                                if rec[j][0] in ("NA","OA"):
                                    orient.extend(quality(li,xyz[h][2],rj) for h in lh[i])
                                if lig[i][0] in ("NA","OA"):
                                    orient.extend(quality(rj,rec[h][2],li) for h in rh[j])
                                if not orient:
                                    continue
                                radius = {"N":1.8,"O":1.7}
                                overlap = radius[lig[i][0][0]] + radius[rec[j][0][0]] - np.linalg.norm(li-rj)
                                if overlap <= 0: continue
                                energy = -.587439 * min(overlap/.7,1)
                                radial += energy
                                penalty -= energy * (1-max(orient))
                                native_contact += 1
                        assert native_contact == int(value["contacts"])
                        max_radial = max(max_radial,abs(radial-float(value["raw_radial"])))
                        max_penalty = max(max_penalty,abs(penalty-float(value["raw_penalty"])))
                        count += 1
    assert max(max_radial,max_penalty) < 1e-9 and count == len(table)
    return dict(poses=count, max_radial_difference=max_radial, max_angular_penalty_difference=max_penalty)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    print(json.dumps(audit(a.output),indent=2))
