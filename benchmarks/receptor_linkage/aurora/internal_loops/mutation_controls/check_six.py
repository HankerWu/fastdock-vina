"""Read the six concentration cells directly from Davis et al. Supplementary Tables 1/4."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import xlrd

HERE = Path(__file__).resolve().parent
SOURCES = HERE.parents[3] / "kinase_state_affinity/source/sources.json"
TERMS = [
    ("AURKA", "VX-680/MK-0457", "Sunitinib"),
    ("CHEK1", "CEP-701", "VX-680/MK-0457"),
    ("S6K1", "Sunitinib", "CEP-701"),
]


def main(path):
    sources = json.loads(SOURCES.read_text())
    books = []
    for name in ("table1.xls", "table4.xls"):
        p = path / name
        assert hashlib.sha256(p.read_bytes()).hexdigest() == sources[name]["sha256"]
        books.append(xlrd.open_workbook(p).sheet_by_index(0))
    table1, table4 = books
    column = table4.row_values(0)
    report = []
    rt = 0.00198720425864083 * 298.15
    for kinase, a, b in TERMS:
        (assay,) = (
            table1.row_values(i)
            for i in range(table1.nrows)
            if table1.cell_value(i, 2) == kinase
        )
        (datum,) = (
            table4.row_values(i)
            for i in range(table4.nrows)
            if table4.cell_value(i, 2) == kinase
        )
        assert assay[:3] == datum[:3] and assay[3] == "NO"
        k1, k2 = (float(datum[column.index(n)]) for n in (a, b))
        report.append(
            dict(
                target=kinase,
                uniprot_or_refseq_accession=assay[0],
                gene=assay[1],
                mutant=assay[3],
                numerator=a,
                denominator=b,
                Kd_nM=[k1, k2],
                RT_log_ratio=rt * math.log(k1 / k2),
            )
        )
    print(
        json.dumps(
            dict(
                temperature_convention_K=298.15,
                RT_kcal_mol=rt,
                source_cells=report,
                closed_sum_kcal_mol=sum(r["RT_log_ratio"] for r in report),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    main(p.parse_args().directory)
