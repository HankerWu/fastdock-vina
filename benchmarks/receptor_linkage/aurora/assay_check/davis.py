"""Check the three frozen ligands against the original Davis Aurora A assay."""

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


def normalized(name):
    return re.sub(r"[^a-z0-9]", "", name.lower())


def main(source, manifest, output):
    settings = json.loads(manifest.read_text())
    for name in ("table1.xls", "table4.xls"):
        digest = hashlib.sha256((source / name).read_bytes()).hexdigest()
        if digest != settings[name]["sha256"]:
            raise ValueError(f"Wrong original workbook: {name}")
    sheet = pd.read_excel(source / "table4.xls")
    rows = sheet[(sheet["Entrez Gene Symbol"] == "AURKA") & (sheet.Kinase == "AURKA")]
    [entry] = rows.to_dict("records")
    if entry["Accession Number"] != "NP_003591.2":
        raise ValueError("Unexpected assay protein")
    aliases = {
        "danusertib": {"danusertib", "pha739358"},
        "at9283": {"at9283"},
        "imatinib": {"imatinib", "sti571", "gleevec", "glivec"},
    }
    with output.open("w") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            ["kinase", "accession", "compound", "drug_column", "status", "kd_nM"]
        )
        for name, names in aliases.items():
            hits = [col for col in entry if normalized(col) in names]
            if len(hits) > 1:
                raise ValueError("Ambiguous drug column")
            col = hits[0] if hits else ""
            value = entry.get(col)
            status = (
                "column_absent"
                if not col
                else "blank"
                if pd.isna(value)
                else "quantified"
            )
            writer.writerow(
                [
                    "AURKA",
                    entry["Accession Number"],
                    name,
                    col,
                    status,
                    value if status == "quantified" else "",
                ]
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    main(args.source, args.manifest, args.output)
