"""Recreate shared docking inputs from the original Kd tables and PDB structures."""

import argparse
import csv
import hashlib
import json
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

import numpy as np
import xlrd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

HERE = Path(__file__).resolve().parent
KINASE = HERE.parents[1] / "kinase_state_affinity"
sys.path.insert(0, str(KINASE / "source"))
from freeze import COMPOUNDS
from prepare import ligand, receptor

STRUCTURES = {"EGFR": ("1M17", "AQ4"), "p38-alpha": ("1A9U", "SB2")}
EXTRA = (
    "AST-487",
    "BIBW-2992",
    "BIRB-796",
    "EXEL-2880/GSK-1363089",
    "Nilotinib",
    "PD-173955",
    "SB-203580",
    "VX-745",
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table_row(sheet, target):
    (row,) = (
        sheet.row_values(i)
        for i in range(1, sheet.nrows)
        if sheet.cell_value(i, 2) == target
    )
    return row


def xyz(lines):
    return np.array([[float(s[i : i + 8]) for i in (30, 38, 46)] for s in lines])


def pocket_audit(work, pdb, component):
    lines = (work / "raw" / (pdb + ".pdb")).read_text().splitlines()
    atoms = [
        s
        for s in lines
        if s.startswith(("ATOM  ", "HETATM")) and s[76:78].strip() not in ("H", "D")
    ]
    ref = xyz(
        [s for s in atoms if s[17:20] == component and s[21] == "A" and s[16] in " A"]
    )
    selected = {
        s for s in atoms if s.startswith("ATOM  ") and s[21] == "A" and s[16] in " A"
    }
    keys = {(s[12:16].strip(), s[22:27].strip()) for s in selected}
    omitted = [
        s for s in atoms if s not in selected and s[17:20] not in (component, "HOH")
    ]
    rebuilt = [
        s
        for s in (work / "prepared" / (pdb + "-protein.pdb")).read_text().splitlines()
        if s.startswith("ATOM")
        and s[76:78].strip() not in ("H", "D")
        and (s[12:16].strip(), s[22:27].strip()) not in keys
    ]

    def near(rows):
        distances = [
            (s[12:27].strip(), float(np.linalg.norm(ref - xyz([s])[0], axis=1).min()))
            for s in rows
        ]
        return [(label, round(d, 3)) for label, d in distances if d < 8]

    return {
        "rebuilt_within_8A": near(rebuilt),
        "omitted_nonwater_within_8A": near(omitted),
    }


def main(work):
    raw, pre = work / "raw", work / "prepared"
    raw.mkdir(parents=True, exist_ok=True)
    pre.mkdir(exist_ok=True)
    for name, item in json.loads((HERE / "sources.json").read_text()).items():
        path = raw / name
        if not path.exists():
            path.write_bytes(urllib.request.urlopen(item["url"], timeout=60).read())
        if digest(path) != item["sha256"]:
            raise ValueError("Source digest mismatch: " + name)
    begin = time.perf_counter()
    with (KINASE / "panel.csv").open() as stream:
        old = {r["compound"]: r for r in csv.DictReader(stream)}
    names = dict(COMPOUNDS)
    graphs = json.loads((raw / "ligands_iso.txt").read_text())
    for name in EXTRA:
        alias = name.split("/")[-1]
        (prop,) = json.loads((raw / (alias + ".json")).read_text())["PropertyTable"][
            "Properties"
        ]
        cid = str(prop["CID"])
        if Chem.MolToInchiKey(Chem.MolFromSmiles(graphs[cid])) != prop["InChIKey"]:
            raise ValueError("Ambiguous compound identity: " + name)
        names[name] = int(cid)
    with tarfile.open(KINASE / "docking-evidence.tar.xz") as archive:
        for row in old.values():
            name = row["cid"] + ".pdbqt"
            (pre / name).write_bytes(archive.extractfile("prepared/" + name).read())
    table1, table4 = [
        xlrd.open_workbook(raw / f"table{n}.xls").sheet_by_index(0) for n in (1, 4)
    ]
    rows, metadata = [], {}
    for target, (pdb, component) in STRUCTURES.items():
        assay = table_row(table4, target)
        ident = table_row(table1, target)
        assert ident[:3] == assay[:3] and ident[3] == "NO"
        t = time.perf_counter()
        receptor(work, pdb, component)
        times = {"receptor_s": time.perf_counter() - t}
        metadata[pdb] = {"times": times, "pocket": pocket_audit(work, pdb, component)}
        active = {
            name: val
            for name, val in zip(table4.row_values(0)[3:], assay[3:])
            if val != ""
        }
        if target == "p38-alpha" and not set(active) <= names.keys():
            raise ValueError("Unmapped measured p38 compound")
        for name in names if target == "p38-alpha" else COMPOUNDS:
            if name not in active:
                continue
            cid = str(names[name])
            smi = graphs[cid]
            if name in old and old[name]["smiles"] != smi:
                raise ValueError("Historical graph mismatch")
            if name not in old:
                t = time.perf_counter()
                ligand(work, {"cid": cid, "smiles": smi})
                times[cid + "_ligand_s"] = time.perf_counter() - t
            mol = max(
                Chem.GetMolFrags(Chem.MolFromSmiles(smi), asMols=True),
                key=lambda m: m.GetNumHeavyAtoms(),
            )
            rows.append(
                {
                    "target": target,
                    "pdb": pdb,
                    "cid": cid,
                    "compound": name,
                    "kd_nm": float(active[name]),
                    "assay_accession": assay[0],
                    "old_panel": int(name in COMPOUNDS),
                    "cocrystal": int(
                        (pdb, name) in {("1M17", "Erlotinib"), ("1A9U", "SB-203580")}
                    ),
                    "heavy": mol.GetNumHeavyAtoms(),
                    "charge": Chem.GetFormalCharge(mol),
                    "scaffold": MurckoScaffold.MurckoScaffoldSmiles(mol=mol),
                    "smiles": Chem.MolToSmiles(mol),
                    "ligand_sha256": digest(pre / (cid + ".pdbqt")),
                }
            )
        metadata[pdb]["receptor_sha256"] = digest(pre / (pdb + ".pdbqt"))
        metadata[pdb]["box_sha256"] = digest(pre / (pdb + ".box"))
    with (work / "panel.csv").open("w") as stream:
        out = csv.DictWriter(stream, list(rows[0]), lineterminator="\n")
        out.writeheader()
        out.writerows(rows)
    metadata["shared_prepare_s_excluding_download"] = time.perf_counter() - begin
    (work / "preparation.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    main(parser.parse_args().work)
