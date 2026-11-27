"""Independent 4LMN redocking with deposited ANP and Mg retained in the receptor."""

import argparse
import json
import csv
import re
import subprocess
import tempfile
import itertools
import sys
from pathlib import Path

import gemmi
import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "pair_gradient/pose"))
import paired
from prepare import base
from sources import crystal_contacts


def rigid_anp(ccd, xyz):
    # Model the phosphate chain as ANP(4-); the crystal constrains heavy atoms.
    block = gemmi.cif.read_file(str(ccd)).sole_block()
    atoms = block.get_mmcif_category("_chem_comp_atom.")
    bonds = block.get_mmcif_category("_chem_comp_bond.")
    ionized = {"O2A", "O2B", "O2G", "O3G"}
    graph, ix = Chem.RWMol(), {}
    for name, el, charge in zip(
        atoms["atom_id"], atoms["type_symbol"], atoms["charge"]
    ):
        if el == "H":
            continue
        a = Chem.Atom(el.title())
        a.SetFormalCharge(-1 if name in ionized else int(charge))
        ix[name] = graph.AddAtom(a)
    if set(ix) != set(xyz):
        raise ValueError("Incomplete measured nucleotide")
    for a, b, order in zip(
        bonds["atom_id_1"], bonds["atom_id_2"], bonds["value_order"]
    ):
        if a in ix and b in ix:
            graph.AddBond(
                ix[a],
                ix[b],
                {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE}[order],
            )
    m = graph.GetMol()
    Chem.SanitizeMol(m)
    assert Chem.GetFormalCharge(m) == -4
    conf = Chem.Conformer(len(ix))
    for name, i in ix.items():
        conf.SetAtomPosition(i, xyz[name])
    m.AddConformer(conf)
    m = Chem.AddHs(m, addCoords=True)
    (setup,) = MoleculePreparation().prepare(m)
    block, ok, error = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(error)
    lines = [
        line for line in block.splitlines() if line.startswith(("ATOM  ", "HETATM"))
    ]
    if sum(line.split()[-1] not in {"HD", "H"} for line in lines) != len(ix):
        raise ValueError("Nucleotide atoms lost during typing")
    return [
        "HETATM" + f"{9100 + i:5d}" + line[11:17] + "ANP A 501" + line[26:]
        for i, line in enumerate(lines)
    ]


def add_crystal_copies(path, structure, center):
    # Include entire rigid units that can contribute to a 24 A box at an
    # 8 A receptor cutoff. Enumeration bounds come from the lattice and the
    # box only; no query-ligand atoms enter the receptor copy.
    lines = [
        l for l in path.read_text().splitlines() if l.startswith(("ATOM  ", "HETATM"))
    ]
    xyz = np.array([base.coords(l) for l in lines])
    heavy = np.array([l.split()[-1] not in ("H", "HD") for l in lines])
    orth = np.array(structure.cell.orth.mat)
    frac = np.array(structure.cell.frac.mat)
    fxyz = xyz @ frac.T
    corners = np.array(
        [
            np.array(center) + offset
            for offset in itertools.product((-20.0, 20.0), repeat=3)
        ]
    )
    fc = corners @ frac.T
    out, copies = [], []
    for op in structure.find_spacegroup().operations():
        rot, tran = np.array(op.rot) / op.DEN, np.array(op.tran) / op.DEN
        image = fxyz @ rot.T + tran
        low = np.ceil(fc.min(axis=0) - image[heavy].max(axis=0)).astype(int)
        high = np.floor(fc.max(axis=0) - image[heavy].min(axis=0)).astype(int)
        for offset in itertools.product(
            *(range(lo, hi + 1) for lo, hi in zip(low, high))
        ):
            current = (image + offset) @ orth.T
            if not np.any(np.all(np.abs(current[heavy] - center) <= 20, axis=1)):
                continue
            copies.append([op.triplet(), offset])
            for line, point in zip(lines, current):
                serial = len(out) + 1
                if serial > 99999:
                    raise ValueError("PDBQT serial overflow")
                out.append(
                    line[:6]
                    + f"{serial:5d}"
                    + line[11:30]
                    + "".join(f"{v:8.3f}" for v in point)
                    + line[54:]
                )
    path.write_text("\n".join(out) + "\n")
    return copies


def prepare(work):
    raw, out = work / "raw", work / "4LMN"
    raw.mkdir(parents=True, exist_ok=True)
    out.mkdir(exist_ok=True)
    source = {}
    for name, url in [
        ("4LMN.pdb", "download/4LMN.pdb"),
        ("EUI.cif", "ligands/download/EUI.cif"),
        ("ANP.cif", "ligands/download/ANP.cif"),
    ]:
        source[name] = base.fetch(raw / name, "https://files.rcsb.org/" + url)
    lines = (raw / "4LMN.pdb").read_text().splitlines()
    observed = base.dep_records(lines, "EUI")
    ligand = np.asarray(list(observed[0].values()))
    meta = base.receptor(lines, out, "4LMN", ligand)
    meta.update(base.ligand(raw / "EUI.cif", observed, out / "ligand.pdbqt"))
    anp = base.dep_records(lines, "ANP")[0]
    cofactor = rigid_anp(raw / "ANP.cif", anp)
    with (out / "receptor.pdbqt").open("a") as io:
        io.write("\n".join(cofactor) + "\n")
    meta.update(
        target="MAP2K1",
        pdb="4LMN",
        chain="A",
        center=ligand.mean(axis=0).tolist(),
        source_sha256=source,
    )
    meta["ca"] = {
        line[22:27].strip() + ":" + line[17:20]: base.coords(line)
        for line in lines
        if line.startswith("ATOM  ")
        and line[21] == "A"
        and line[16] in " A"
        and line[12:16].strip() == "CA"
    }
    meta["crystal_contacts"] = crystal_contacts(
        gemmi.read_structure(str(raw / "4LMN.pdb")), "A", "EUI"
    )
    meta["nearest_ANP_heavy_A"] = min(
        float(np.linalg.norm(ligand - v, axis=1).min()) for v in anp.values()
    )
    meta["crystal_copies"] = add_crystal_copies(
        out / "receptor.pdbqt",
        gemmi.read_structure(str(raw / "4LMN.pdb")),
        meta["center"],
    )
    if (
        meta["metals"] != ["MG"]
        or meta["missing_atoms"]
        or len(meta["crystal_copies"]) < 2
    ):
        raise ValueError(
            f"4LMN input gate: {meta['metals']}, {meta['missing_atoms']}, {meta['crystal_contacts']}"
        )
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")


def fixed_pool(work, meta, parent, fixed):
    rows = []
    cells = list(csv.DictReader((work / "cells.csv").open()))
    with tempfile.TemporaryDirectory() as temp:
        pose = Path(temp) / "pose.pdbqt"
        for cell in cells:
            if cell["arm"] != "parent":
                continue
            pool = (work / "runs" / (cell["key"] + ".pdbqt")).read_text()
            poses = re.findall(r"^MODEL .*?^ENDMDL", pool, re.M | re.S)
            if len(poses) != int(cell["modes"]):
                raise ValueError("Incomplete parent pool")
            for rank, block in enumerate(poses, 1):
                pose.write_text(block + "\n")
                for reference in (
                    ("common_zero", "self") if rank == 1 else ("common_zero",)
                ):
                    cmd = [
                        "--receptor",
                        work / "4LMN/receptor.pdbqt",
                        "--ligand",
                        work / "4LMN/ligand.pdbqt",
                        "--pose_coordinates",
                        pose,
                        "--score_only",
                        "--cpu",
                        "1",
                        "--seed",
                        "817139",
                    ]
                    if reference == "common_zero":
                        cmd.extend(["--unbound_energy", "0"])
                    for axis, x in zip("xyz", meta["4LMN"]["center"]):
                        cmd.extend(["--center_" + axis, str(x), "--size_" + axis, "24"])
                    for arm, binary in [("parent", parent), ("fixed", fixed)]:
                        out = subprocess.check_output([binary, *cmd], text=True)
                        values = re.findall(r":\s*([-+.\deE]+)\s+\(kcal/mol\)", out)
                        if len(values) != 10:
                            raise ValueError("Incomplete energy components")
                        rows.append(
                            dict(
                                cap=cell["cap"],
                                parent_rank=rank,
                                reference=reference,
                                arm=arm,
                                total=values[0],
                                inter=values[1],
                                intra=values[4],
                                unbound=values[-1],
                            )
                        )
    paired.write(work / "fixed_pools.csv", rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=["prepare", "run", "fixed-pool"])
    p.add_argument("work", type=Path)
    for arm in ("upstream", "counts", "parent", "fixed"):
        p.add_argument("--" + arm, type=Path)
    args = p.parse_args()
    work = args.work.resolve()
    if args.action == "prepare":
        return prepare(work)
    meta = {"4LMN": json.loads((work / "4LMN/metadata.json").read_text())}
    if args.action == "fixed-pool":
        return fixed_pool(work, meta, args.parent, args.fixed)
    (work / "runs").mkdir(exist_ok=True)
    paired.RADII["Mg"] = 1.73
    runs, modes = [], []
    for cap, order in [
        (3000, ("fixed", "upstream", "parent")),
        (12000, ("parent", "upstream", "fixed")),
    ]:
        for arm in order:
            r, ps = paired.execute(
                work,
                meta,
                "4LMN",
                "4LMN",
                817139,
                cap,
                arm,
                getattr(args, arm),
                False,
                args.counts,
            )
            runs.append(r)
            modes.extend(ps)
    paired.write(work / "cells.csv", runs)
    paired.write(work / "modes.csv", modes)


if __name__ == "__main__":
    main()
