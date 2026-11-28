"""Independent raw-PDBQT, primary-assay and tie-correct rank verification."""

import argparse
import csv
from collections import defaultdict
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from math import isclose, sqrt
from pathlib import Path
from statistics import correlation, median
import xlrd

SOURCE = Path(__file__).resolve().parents[2] / "geometry_comparison/affinity"


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def original_labels(raw, panel):
    data = xlrd.open_workbook(str(raw / "table4.xls")).sheet_by_index(0)
    ids = xlrd.open_workbook(str(raw / "table1.xls")).sheet_by_index(0)
    kd = {data.cell_value(i, 2): data.row_values(i) for i in range(1, data.nrows)}
    wt = {ids.cell_value(i, 2): ids.row_values(i) for i in range(1, ids.nrows)}
    assay_records = defaultdict(set)
    for path in raw.glob("*_Ki.json"):
        for a in json.loads(path.read_text())["activities"]:
            if (
                a["standard_relation"],
                a["standard_units"],
                a["standard_type"],
                a["target_organism"],
            ) == ("=", "nM", "Ki", "Homo sapiens"):
                assay_records[a["assay_chembl_id"], a["molecule_chembl_id"]].add(
                    Decimal(a["standard_value"])
                )
    for r in panel:
        if r["endpoint"] == "Kd":
            target = r["assay"][6:]
            assert wt[target][3] == "NO" and wt[target][0] == kd[target][0]
            value = Decimal(str(kd[target][data.row_values(0).index(r["name"])]))
        else:
            (value,) = assay_records[r["assay"], r["compound"]]
        assert value == Decimal(r["value_nM"])


def rank(values):
    return [sum(x < v for x in values) + (1 + values.count(v)) / 2 for v in values]


def sign(v):
    return int(v > 0) - int(v < 0)


def independent(work, results):
    panel = read_csv(SOURCE / "panel.csv")
    original_labels(work / "raw", panel)
    sources = json.loads((SOURCE / "sources.json").read_text())
    for file, source in sources.items():
        assert hashlib.sha256((work / "raw" / file).read_bytes()).hexdigest() == source["sha256"]
    measured = {(r["assay"], r["compound"]): Decimal(r["value_nM"]) for r in panel}
    records = read_csv(results / "runs.csv")
    predictions = defaultdict(list)
    raw_outputs = set()
    for r in records:
        p = (
            work
            / "outputs"
            / r["engine"]
            / r["storage_stage"]
            / f"{r['receptor']}_{r['compound']}_{r['seed']}.pdbqt"
        )
        key = r["stage"], r["engine"], r["assay"], r["receptor"], r["compound"]
        valid = r["valid"] == "True"
        if valid:
            blob = p.read_bytes()
            assert hashlib.sha256(blob).hexdigest() == r["output_sha256"]
            text = blob.decode("ascii").splitlines()
            assert text[0] == "MODEL 1" and text[-1] == "ENDMDL"
            assert sum(line.startswith("MODEL ") for line in text) == int(r["modes"])
            score = Decimal(text[1].split(":")[1].split()[0])
            assert score == Decimal(r["top_score"])
            raw_outputs.add(p)
            predictions[key].append(score)
        else:
            predictions[key].append(None)
    scores = {k: median(v) for k, v in predictions.items() if all(x is not None for x in v)}
    summary = read_csv(results / "assays.csv")
    pairs = read_csv(results / "pairs.csv")
    for row in pairs:
        prefix = tuple(row[k] for k in ("stage", "engine", "assay", "receptor"))
        a, b = row["stronger"], row["weaker"]
        assert measured[row["assay"], a] < measured[row["assay"], b]
        if (*prefix, a) in scores and (*prefix, b) in scores:
            diff = scores[*prefix, b] - scores[*prefix, a]
            assert diff == Decimal(row["weaker_minus_stronger"])
            assert Fraction(1 + sign(diff), 2) == Fraction(row["credit"])
    for s in summary:
        prefix = tuple(s[k] for k in ("stage", "engine", "assay", "receptor"))
        vals = [(measured[s["assay"], k[-1]], v) for k, v in scores.items() if k[:-1] == prefix]
        assert int(s["observed_compounds"]) == len(vals)
        if not vals:
            continue
        x, y = map(list, zip(*vals))
        xy = [
            (sign(a - c), sign(b - d), max(a, c) / min(a, c)) for (a, b), (c, d) in combinations(vals, 2)
        ]
        if len(set(y)) > 1 and len(set(x)) > 1:
            tau = sum(u * v for u, v, _ in xy) / sqrt(
                sum(u != 0 for u, _, _ in xy) * sum(v != 0 for _, v, _ in xy)
            )
            assert isclose(float(s["tau_b"]), tau, abs_tol=1e-10)
            assert isclose(float(s["rho"]), correlation(rank(x), rank(y)), abs_tol=1e-10)
        for t in (1, 3, 10):
            eligible = [u * v for u, v, fold in xy if u and fold >= t]
            assert len(eligible) == int(s[f"pairs_{t}x"])
            assert sum((1 + d) / 2 for d in eligible) == float(s[f"correct_{t}x"])
    return dict(
        source_hashes=len(sources),
        primary_assay_labels=len(panel),
        unique_raw_top_scores=len(raw_outputs),
        planned_runs=len(records),
        rank_summaries=len(summary),
        verified_pair_rows=len(pairs),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    print(json.dumps(independent(args.work, args.results), indent=2))
