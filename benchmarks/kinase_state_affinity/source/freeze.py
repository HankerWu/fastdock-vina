"""Extract Davis Suppl. Table 4 exact Kd from three unmodified kinase assay rows."""

import csv
import json
import sys
from pathlib import Path


COMPOUNDS = {
    "AT-7519": 11338033,
    "BMS-387032/SNS-032": 3025986,
    "CEP-701": 126565,
    "Flavopiridol": 5287969,
    "LY-317615": 176167,
    "Sorafenib": 216239,
    "Staurosporine": 44259,
    "Dasatinib": 3062316,
    "Erlotinib": 176870,
    "Gefitinib": 123631,
    "SKI-606": 5328940,
    "Sunitinib": 5329102,
    "R406": 11984591,
}
RECEPTORS = {
    "CDK2": [("1HCK", "ATP"), ("1G5S", "I17")],
    "GSK3B": [("1Q3D", "STU"), ("1Q5K", "TMU")],
    "SRC": [("1Y57", "MPZ"), ("2SRC", "ANP")],
}
SEEDS = (271828, 314159)


def freeze(work):
    import xlrd
    from rdkit import Chem

    sheet = xlrd.open_workbook(work / "raw/table4.xls").sheet_by_index(0)
    headings = sheet.row_values(0)
    structures = json.loads((work / "raw/ligands_iso.txt").read_text())
    keys = {
        d["CID"]: d["InChIKey"]
        for d in json.loads((work / "raw/identities.json").read_text())[
            "PropertyTable"
        ]["Properties"]
    }
    with (work / "panel.csv").open("w") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["target", "accession", "compound", "cid", "kd_nm", "smiles", "inchi_key"]
        )
        for target in RECEPTORS:
            (row,) = [
                sheet.row_values(i)
                for i in range(1, sheet.nrows)
                if sheet.cell_value(i, 2) == target
            ]
            for name, cid in COMPOUNDS.items():
                value = row[headings.index(name)]
                if value == "":
                    continue
                value = float(value)
                assert value > 0
                smiles = structures[str(cid)]
                key = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
                if keys[cid] != key:
                    raise ValueError(f"{cid}: {keys[cid]} != {key}")
                writer.writerow([target, row[0], name, cid, value, smiles, key])


if __name__ == "__main__":
    freeze(Path(sys.argv[1]))
