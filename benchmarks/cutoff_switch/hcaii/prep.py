"""Prepare hCA II ligands from CCD graph, protein from experimental ATOM records."""

import argparse
import concurrent.futures
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import gemmi
import numpy as np
from download import fetch_sources
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

CASES = {
    "6WQ4": ("U7A", "34"),
    "6WQ5": ("U7G", "41"),
    "6WQ7": ("U6V", "42"),
    "6WQ8": ("U77", "46"),
    "6WQ9": ("U7J", "48"),
}
RECS = ["6WQ4", "6WQ8"]


def records(raw, p):
    return (raw / (p + ".pdb")).read_text().splitlines()


def xyz(l):
    return np.array([float(l[30 + i : 38 + i]) for i in (0, 8, 16)])


def receptor(work, p, box_size):
    P = work / "prepared"
    lines = records(work / "raw", p)
    protein = [
        l
        for l in lines
        if l.startswith("ATOM  ")
        and l[21] == "A"
        and l[16] in " A"
        and l[76:78].strip() not in ["H", "D"]
    ]
    zinc = next(
        l
        for l in lines
        if l.startswith("HETATM") and l[17:20].strip() == "ZN" and l[21] == "A"
    )
    clean = P / f"{p}_protein.pdb"
    clean.write_text("\n".join(l[:16] + " " + l[17:] for l in protein) + "\nEND\n")
    c = xyz(zinc)
    donors = np.array(
        [
            xyz(
                next(l for l in protein if int(l[22:26]) == n and l[12:16].strip() == a)
            )
            for n, a in [(94, "NE2"), (96, "NE2"), (119, "ND1")]
        ]
    )
    vec = c - donors
    vec /= np.linalg.norm(vec, axis=1)[:, None]
    axis = vec.sum(axis=0)
    center = c + 5 * axis / np.linalg.norm(axis)
    (P / f"{p}_box.txt").write_text(
        "".join(f"center_{a} = {v:.3f}\n" for a, v in zip("xyz", center))
        + "".join(f"size_{a} = {box_size:g}\n" for a in "xyz")
    )
    fn = P / f"{p}_receptor.pdbqt"
    # Histidine tautomer assignments leave the observed Zn-coordinating nitrogen unprotonated.
    cmd = [
        sys.executable,
        "-m",
        "meeko.cli.mk_prepare_receptor",
        "--read_pdb",
        str(clean),
        "--set_template",
        "A:94,96=HID,A:119=HIE",
        "-p",
        str(fn),
    ]
    o = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (P / f"{p}_receptor.log").write_text(o.stdout)
    if o.returncode:
        raise RuntimeError(str(cmd) + o.stdout[-2000:])
    qt = fn.read_text()
    serial = 1 + max(
        int(l[6:11]) for l in qt.splitlines() if l.startswith(("ATOM", "HETATM"))
    )
    zn = zinc[:6] + f"{serial:5d}" + zinc[11:54] + "  1.00  0.00    +2.000 Zn"
    fn.write_text(qt.rstrip() + "\n" + zn + "\n")


def ca(raw, p):
    return {
        int(l[22:26]): xyz(l)
        for l in records(raw, p)
        if l.startswith("ATOM  ")
        and l[12:16].strip() == "CA"
        and l[21] == "A"
        and l[16] in " A"
    }


def transform(raw, src, dst, X):
    a, b = ca(raw, src), ca(raw, dst)
    ids = sorted(a.keys() & b.keys())
    A = np.array([a[i] for i in ids])
    B = np.array([b[i] for i in ids])
    aa = A.mean(0)
    bb = B.mean(0)
    u, s, vt = np.linalg.svd((A - aa).T @ (B - bb))
    u[:, -1] *= np.sign(np.linalg.det(u @ vt))
    rotation = u @ vt
    return (
        (X - aa) @ rotation + bb,
        len(ids),
        np.sqrt(np.mean(np.sum(((A - aa) @ rotation + bb - B) ** 2, axis=1))),
    )


def ligand(work, p):
    R = work / "raw"
    P = work / "prepared"
    c, num = CASES[p]
    b = gemmi.cif.read_file(str(R / f"{c}.cif")).sole_block()
    atoms = b.find("_chem_comp_atom.", ["atom_id", "type_symbol", "charge"])
    m = Chem.RWMol()
    idx = {}
    names = []
    for name, el, chg in atoms:
        if el == "H":
            continue
        atom = Chem.Atom(el.capitalize())
        atom.SetFormalCharge(int(chg))
        idx[name] = m.AddAtom(atom)
        names.append(name)
    for a, d, t in b.find(
        "_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]
    ):
        if a in idx and d in idx:
            m.AddBond(
                idx[a],
                idx[d],
                {
                    "SING": Chem.BondType.SINGLE,
                    "DOUB": Chem.BondType.DOUBLE,
                    "TRIP": Chem.BondType.TRIPLE,
                    "AROM": Chem.BondType.AROMATIC,
                }[t],
            )
    m = m.GetMol()
    Chem.SanitizeMol(m)
    # Fixed pH 7.5 model: zinc-bound sulfonamide anion, free aliphatic amines protonated and acids ionized.
    for pat, q in [
        ("[N;H2;+0:1]S(=O)(=O)", -1),
        ("[N;X3;+0;!$(NS);!$(NC=O);!$(N#*)]", 1),
        ("[O;H1;+0]C=O", -1),
    ]:
        for match in m.GetSubstructMatches(Chem.MolFromSmarts(pat)):
            atom = m.GetAtomWithIdx(match[0])
            atom.SetFormalCharge(q)
            atom.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(m)
    # Generate docking input from independently embedded conformer; reference coordinates are read separately.
    h = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = 9033
    if AllChem.EmbedMolecule(h, params) != 0:
        raise ValueError(p)
    AllChem.MMFFOptimizeMolecule(h, maxIters=1000)
    setup = MoleculePreparation().prepare(h)[0]
    qt, ok, reason = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(reason)
    (P / f"{p}_ligand.pdbqt").write_text(qt)
    smi = Chem.MolToSmiles(m)
    print(p, smi, flush=True)
    serial = {}
    for l in qt.splitlines():
        if l.startswith("REMARK INDEX MAP"):
            aa = list(map(int, l.split()[3:]))
            serial.update({a - 1: d for a, d in zip(aa[::2], aa[1::2])})
    graph_sdf = P / f"{p}_graph.sdf"
    with Chem.SDWriter(str(graph_sdf)) as f:
        f.write(Chem.RemoveHs(h))
    crystal = {
        l[12:16].strip(): xyz(l)
        for l in records(R, p)
        if l.startswith("HETATM")
        and l[17:20] == c
        and l[21] == "A"
        and int(l[22:26]) == 303
        and l[16] in " A"
        and l[76:78].strip() not in ["H", "D"]
    }
    assert set(crystal) == set(names), (p, set(crystal) ^ set(names))
    labels = ET.parse(R / "paper.xml").getroot().find(".//table-wrap[@id='tbl1']")
    rows = [["".join(x.itertext()) for x in tr] for tr in labels.findall(".//tr")]
    assert next(row for row in rows if row[0] == "cmpd")[6] == "CA II"
    ki = float(next(row for row in rows if row[0] == num)[6])
    ref = {}
    for rec in dict.fromkeys([*RECS, p]):
        coords, nca, rmsd = transform(R, p, rec, np.array([crystal[n] for n in names]))
        ref[rec] = {
            "coords": np.round(coords, 6).tolist(),
            "aligned_CA_count": nca,
            "aligned_CA_rmsd": float(rmsd),
        }
    return p, {
        "ccd": c,
        "compound": num,
        "Ki_nM": ki,
        "smiles": smi,
        "names": names,
        "qt_serial": [serial[i] for i in range(len(names))],
        "references": ref,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--box-size", type=float, default=36.0)
    args = parser.parse_args()
    work = args.work.resolve()
    fetch_sources(work / "raw")
    prepared = work / "prepared"
    prepared.mkdir(parents=True, exist_ok=False)
    with concurrent.futures.ThreadPoolExecutor(3) as pool:
        receptors = [pool.submit(receptor, work, pdb, args.box_size) for pdb in CASES]
        ligands = dict(ligand(work, pdb) for pdb in CASES)
        for result in receptors:
            result.result()
    for ligand_pdb, ligand_data in ligands.items():
        for receptor_pdb, ref in ligand_data["references"].items():
            box = [
                float(line.split("=")[1])
                for line in (prepared / f"{receptor_pdb}_box.txt")
                .read_text()
                .splitlines()
            ]
            margin = np.min(
                np.array(box[3:]) / 2 - np.abs(np.array(ref["coords"]) - box[:3])
            )
            ref["minimum_box_margin"] = round(float(margin), 5)
            if margin < 0:
                print(
                    f"Reference outside box: {receptor_pdb}/{ligand_pdb}: {margin:.3f} A"
                )
    (prepared / "metadata.json").write_text(json.dumps(ligands, indent=2) + "\n")


if __name__ == "__main__":
    main()
