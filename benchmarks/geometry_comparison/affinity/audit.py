"""Independent raw-output/assay audit using only standard Python and xlrd."""
import argparse
import csv
import itertools
import json
from pathlib import Path

import xlrd


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def audit(work, raw, seeds):
    rows = read(work / "panel.csv")
    raw_ki = {}
    for file in raw.glob("*_Ki.json"):
        for rec in json.loads(file.read_text())["activities"]:
            if rec["standard_relation"] == "=" and rec["standard_units"] == "nM":
                raw_ki[rec["assay_chembl_id"], rec["molecule_chembl_id"]] = float(rec["standard_value"])
    data = xlrd.open_workbook(raw / "table4.xls").sheet_by_index(0)
    davis = {data.cell_value(i, 2): data.row_values(i) for i in range(1, data.nrows)}
    kd_heads = data.row_values(0)
    for row in rows:
        if row["endpoint"] == "Ki":
            value = raw_ki[row["assay"], row["compound"]]
        else:
            value = float(davis[row["assay"].removeprefix("Davis-")][kd_heads.index(row["name"])])
        assert value == float(row["value_nM"])

    pools = work / "pools"
    seen = {}
    for pose in pools.glob("*.pdbqt"):
        state, compound, seed, arm = pose.stem.split("_")
        if int(seed) not in seeds:
            continue
        assert arm in ("armijo", "geometry")
        text = pose.read_text()
        scores = [float(l.split()[3]) for l in text.splitlines() if l.startswith("REMARK VINA RESULT:")]
        assert min(scores) == scores[0]
        seen[state, compound, arm, int(seed)] = scores[0]
    assert len(seen) == len(read(work / "cells.csv"))
    for summary in read(work / "summary.csv"):
        subset = [r for r in rows if r["assay"] == summary["assay"]]
        predicted = {}
        for r in subset:
            values = [seen[summary["receptor"], r["compound"], summary["arm"], seed] for seed in seeds]
            assert len(values) == 2
            predicted[r["compound"]] = sum(values) / 2
        pairs, correct = 0, 0
        for a, b in itertools.permutations(subset, 2):
            if float(a["value_nM"]) >= float(b["value_nM"]):
                continue
            if float(b["value_nM"]) / float(a["value_nM"]) < float(summary["fold_min"]):
                continue
            pairs += 1
            delta = predicted[a["compound"]] - predicted[b["compound"]]
            correct += (delta < 0) + 0.5 * (delta == 0)
        assert (pairs, correct) == (int(summary["pairs"]), float(summary["correct"]))
    return dict(raw_docking_outputs=len(seen), experiment_rows=len(rows), verified_summary_rows=len(read(work / "summary.csv")))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("raw", type=Path)
    p.add_argument("seeds", type=int, nargs=2)
    args = p.parse_args()
    print(json.dumps(audit(args.work, args.raw, args.seeds), indent=2))
