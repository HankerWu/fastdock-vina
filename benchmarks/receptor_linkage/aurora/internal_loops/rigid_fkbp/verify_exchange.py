"""Check the published conditional-score bounds and read the six-cell cycle from Davis Table 4."""

import argparse
from decimal import Decimal as D, getcontext
import json
from pathlib import Path
import subprocess
import sys
import xlrd
from run_panel import ROOT, sha

getcontext().prec = 40


def verify(source):
    bench = ROOT / "benchmarks"
    magn = bench / "affinity_exchange/joint_reference/local_integral/magnitude"
    evidence = subprocess.check_output(
        [sys.executable, str(magn / "verify.py"), str(magn / "results")], text=True
    )
    manifest = json.loads(
        (bench / "kinase_state_affinity/source/sources.json").read_text()
    )["table4.xls"]
    assert sha(source) == manifest["sha256"]
    sheet = xlrd.open_workbook(str(source)).sheet_by_index(0)
    header = sheet.row_values(0)
    cycle = (
        ("AURKA", "VX-680/MK-0457", "Sunitinib"),
        ("CHEK1", "CEP-701", "VX-680/MK-0457"),
        ("S6K1", "Sunitinib", "CEP-701"),
    )
    terms = []
    for target, a, b in cycle:
        (row,) = (row for row in sheet.get_rows() if row[2].value == target)
        ka, kb = (D(str(row[header.index(k)].value)) for k in (a, b))
        terms.append(
            dict(
                target=target,
                a=a,
                b=b,
                kd_a_nM=str(ka),
                kd_b_nM=str(kb),
                log_ratio_kcal_mol=str(
                    D(".00198720425864083") * D("298.15") * (ka / kb).ln()
                ),
            )
        )
    return dict(
        source=manifest,
        checks=evidence.strip(),
        terms=terms,
        cycle_kcal_mol=str(sum(D(t["log_ratio_kcal_mol"]) for t in terms)),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("davis_table4_xls", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    args.out.write_text(json.dumps(verify(args.davis_table4_xls), indent=2) + "\n")
