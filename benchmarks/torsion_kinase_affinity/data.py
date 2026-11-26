"""Freeze named, measured exact Kd rows and prepare their independent 3D inputs."""

import csv
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib
import sys
import xlrd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def build(work, davis_raw, repo, pose):
    sys.path.insert(0, str(repo / "benchmarks/kinase_state_affinity/source"))
    import freeze, prepare

    ligands = json.loads((davis_raw / "ligands_iso.txt").read_text())
    identities = {
        r["CID"]: r["InChIKey"]
        for r in json.loads((davis_raw / "identities.json").read_text())[
            "PropertyTable"
        ]["Properties"]
    }
    compounds = {
        n: (cid, ligands[str(cid)], identities[cid])
        for n, cid in freeze.COMPOUNDS.items()
    }
    for name, v in json.loads((work / "names.json").read_text()).items():
        (r,) = v["record"]
        compounds[name] = r["CID"], r["SMILES"], r["InChIKey"]
    files = ["table1.xls", "table4.xls"]
    provenance = json.loads(
        (repo / "benchmarks/kinase_state_affinity/source/sources.json").read_text()
    )
    for file in files:
        assert (
            hashlib.sha256((davis_raw / file).read_bytes()).hexdigest()
            == provenance[file]["sha256"]
        )
    sheet1, sheet4 = [
        xlrd.open_workbook(davis_raw / file).sheet_by_index(0) for file in files
    ]
    panel = []
    for t in ("AURKA", "CHEK1"):
        (assay,) = [r for r in sheet1.get_rows() if r[2].value == t]
        assert assay[3].value == "NO"
        (values,) = [r for r in sheet4.get_rows() if r[2].value == t]
        assert values[0].value == assay[0].value
        for col, name in enumerate(sheet4.row_values(0)):
            if name not in compounds or not values[col].value:
                continue
            cid, smi, key = compounds[name]
            m = Chem.MolFromSmiles(smi)
            if Chem.MolToInchiKey(m) != key:
                raise ValueError((name, cid, key))
            parts = Chem.GetMolFrags(m, asMols=True)
            parent = max(parts, key=lambda p: p.GetNumHeavyAtoms())
            counterions = sorted(Chem.MolToSmiles(x) for x in parts if x is not parent)
            expected = {"R406": ["O=S(=O)(O)c1ccccc1"], "JNJ-28312141": ["Cl"]}.get(
                name, []
            )
            if counterions != expected:
                raise ValueError((name, counterions, expected))
            panel.append(
                dict(
                    target=t,
                    accession=assay[0].value,
                    compound=name,
                    cid=str(cid),
                    kd_nm=float(values[col].value),
                    smiles=Chem.MolToSmiles(parent),
                    source_smiles=smi,
                    inchi_key=key,
                    heavy=parent.GetNumHeavyAtoms(),
                    scaffold=Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(parent)),
                )
            )
    with (work / "panel.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=panel[0])
        w.writeheader()
        w.writerows(panel)
    for s in ("2E9P", "2YWP"):
        p = work / "prepared"
        (p / (s + ".pdbqt")).write_bytes((pose / s / "receptor.pdbqt").read_bytes())
        (p / (s + ".json")).write_text(
            json.dumps(
                {
                    "center": json.loads((pose / s / "preparation.json").read_text())[
                        "center"
                    ]
                }
            )
            + "\n"
        )
    rows = {r["cid"]: r for r in panel}

    def gen(r):
        prepare.ligand(work, r)

    with ThreadPoolExecutor(2) as e:
        list(e.map(gen, rows.values()))
    print("validation", len(panel), len(rows), flush=True)


if __name__ == "__main__":
    build(*(Path(a).resolve() for a in sys.argv[1:]))
