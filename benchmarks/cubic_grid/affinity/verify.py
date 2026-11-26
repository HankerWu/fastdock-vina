"""Independent, direct extraction of Davis Kd, pool energies, rank counts and controls."""

import csv
import hashlib
import json
import re
import statistics
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import xlrd

work, table1, table4 = map(Path, sys.argv[1:])
for path, expected in [
    (table1, "69c2957c1fa1b5ba032721031af2fffe1eb896acd0723829e4500da8f75779c9"),
    (table4, "e47821b2e959ff22aa6d1613b236a63fb530f0da19cfa8d45408598cf27498d7"),
]:
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected

s, t = [xlrd.open_workbook(p).sheet_by_index(0) for p in (table1, table4)]
panel = list(csv.DictReader((work / "panel.csv").open()))
summary = {
    r["target"]: r
    for r in csv.DictReader((work / "summary.csv").open())
    if r["seed"] == "mean"
}
cells = {
    (r["target"], r["cid"], r["arm"], r["seed"]): r
    for r in csv.DictReader((work / "cells.csv").open())
}
all_results = {}
for target, state in [("FAK", "6YT6"), ("SYK", "6ZCY")]:
    q = [r for r in panel if r["target"] == target]
    (i,) = [i for i in range(s.nrows) if s.cell_value(i, 2) == target]
    (h,) = [i for i in range(t.nrows) if t.cell_value(i, 2) == target]
    assert s.cell_value(i, 3) == "NO" and s.cell_value(i, 0) == t.cell_value(h, 0)
    energies = defaultdict(list)
    labels = {}
    for r in q:
        labels[r["cid"]] = float(t.cell_value(h, t.row_values(0).index(r["name"])))
        assert labels[r["cid"]] == float(r["kd_nm"])
        for a in ("native", "cubic"):
            for seed in ("76123", "76231"):
                path = work / "pools" / f"{state}_{r['cid']}_{seed}_{a}.pdbqt"
                pool = path.read_text()
                score = float(
                    re.search(r"REMARK VINA RESULT:\s*([-\d.]+)", pool).group(1)
                )
                cell = cells[target, r["cid"], a, seed]
                assert abs(score - float(cell["score"])) < 1e-8
                energies[r["cid"], a].append(score)
                energies[r["cid"], a + "-readback"].append(float(cell["readback"]))
    pairs = [(u, v) for u in labels for v in labels if labels[v] >= 3 * labels[u]]
    totals = {}
    for arm in ("native", "cubic", "native-readback", "cubic-readback"):
        total = 0
        for u, v in pairs:
            diff = statistics.mean(energies[v, arm]) - statistics.mean(energies[u, arm])
            total += 1 if diff > 1e-8 else 0.5 if abs(diff) <= 1e-8 else 0
        totals[arm] = total
    assert totals["native"] == float(summary[target]["native"]) and totals[
        "cubic"
    ] == float(summary[target]["cubic"])
    assert len(pairs) == int(summary[target]["pairs"])
    all_results[target] = {
        "exact_assay_entries": sum(bool(v) for v in t.row_values(h)[3:]),
        "checked_entries": len(q),
        "unique_pairs": len(pairs),
        "pair_counts": totals,
        "all_cids_distinct": len(q) == len(labels),
    }
(work / "independent.json").write_text(json.dumps(all_results, indent=2) + "\n")
print(json.dumps(all_results))
