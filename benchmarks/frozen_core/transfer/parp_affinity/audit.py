"""Audit PARP1 seeded top scores, pool invariance, source labels and expected pair counts."""

import argparse
import csv
from itertools import permutations
from pathlib import Path
import statistics
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent


def verify(work, parent, table_dir):
    table = ET.parse(HERE / "table1.xml")
    levels = {}
    for row in table.findall(".//tbody/tr"):
        cells = list(row)
        label = "".join(cells[0].itertext()).casefold()
        values = "".join(cells[4].itertext()).split("±")
        levels[label] = tuple(map(float, values))
    manual = {
        "talazoparib": (0.6, 0.1),
        "olaparib": (0.8, 0.1),
        "niraparib": (5.4, 1.4),
        "veliparib": (6.2, 3.5),
    }
    assert levels == manual
    comparisons = [
        (a, b) for a, b in permutations(levels, 2) if levels[b][0] / levels[a][0] >= 3
    ]
    assert len(comparisons) == 4
    data = {}
    unchanged = 0
    for p in (work / "pools").glob("*.pdbqt"):
        compound, rec, seed = p.stem.split("-")
        text = p.read_text()
        model = text.split("ENDMDL")[0].splitlines()
        native = float(
            next(l for l in model if l.startswith("REMARK VINA RESULT:")).split()[3]
        )
        contact = float(
            next(
                l for l in model if l.startswith("REMARK VINA CONTACT ROTOR SCORE:")
            ).split()[5]
        )
        data[compound, rec, seed] = (native, contact)
        old = {"olaparib": "7KK4", "talazoparib": "4UND"}.get(compound)
        if old and old != rec and (parent / "pools").exists():
            reference = parent / "pools" / f"{old}-{rec}-{seed}.pdbqt"
            stream = "".join(
                l
                for l in text.splitlines(keepends=True)
                if not l.startswith("REMARK VINA CONTACT ROTOR SCORE:")
            )
            assert reference.read_text() == stream
            unchanged += 1
    assert len(data) == 16
    if (parent / "pools").exists():
        assert unchanged == 4
    with (table_dir / "summary.csv").open() as stream:
        ref = list(csv.DictReader(stream))
    for state in ("4UND", "7KK4", "mean_receptors"):
        recs = ("4UND", "7KK4") if state == "mean_receptors" else [state]
        for i, method in enumerate(("native", "contact")):
            g = {
                c: statistics.mean(
                    statistics.mean(data[c, r, k][i] for k in ("98207", "98213"))
                    for r in recs
                )
                for c in levels
            }
            observed = sum((g[a] < g[b]) + 0.5 * (g[a] == g[b]) for a, b in comparisons)
            (row,) = [
                z for z in ref if z["receptor"] == state and z["method"] == method
            ]
            assert (int(row["pairs"]), float(row["correct"])) == (
                len(comparisons),
                observed,
            )
    print(
        f"4 assay Kd values, 16 native top records, {unchanged} unannotated pool equivalents independently verified"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("parent", type=Path)
    p.add_argument("--tables", type=Path, default=HERE)
    a = p.parse_args()
    verify(a.work, a.parent, a.tables)
