"""Apply the registered predictive gates to the four real source arms."""
import argparse
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STRUCTURE = HERE.parent / "fork_structure"


def read(file):
    with file.open() as stream:
        return list(csv.DictReader(stream))


def decision():
    sections = {}
    for dataset in ("fixed", "external"):
        entry = {}
        for row in read(STRUCTURE / dataset / "summary.csv"):
            if row["target"] == "ALL":
                assert int(row["scheduled"]) == int(row["valid"])
                entry[row["cohort"], row["arm"]] = row
        sections[dataset] = entry
    tests = {}
    for dataset, cohort, metric, strictly in [
        ("fixed", "reciprocal", "mean_first", True), ("fixed", "reciprocal", "mean_best", True),
        ("fixed", "reciprocal", "p95", False), ("external", "cognate", "mean_first", False),
        ("external", "reciprocal", "mean_first", False), ("external", "reciprocal", "mean_best", False),
    ]:
        section = sections[dataset]
        x, y = [float(section[cohort, arm][metric]) for arm in ("stable","candidate")]
        tests[f"{dataset}/{cohort}/{metric}"] = dict(start=x, candidate=y, pass_gate=y<x if strictly else y<=x)
    arm_pair_credit = {}
    for arm in ("official", "parent", "stable", "candidate"):
        rr = [r for r in read(HERE/"affinity/assays.csv") if r["stage"] == "nominal" and r["engine"] == arm]
        assert len(rr) == 9
        n = sum(int(r["pairs_3x"]) for r in rr)
        assert n == sum(int(r["expected_3x"]) for r in rr) == 187
        arm_pair_credit[arm] = sum(float(r["correct_3x"]) for r in rr)
    old, new = [arm_pair_credit[k] for k in ("stable", "candidate")]
    tests["nominal_assay_3x_credit"] = dict(start=old, candidate=new, pass_gate=new>=old)
    return dict(round_start="0449b1d0bc3adec768dd9a9f243f3ce5a29cef58",
                candidate="1c672a8e2ce12c7a16e20172a3315ea5a19c23d8", gates=tests,
                nominal_assay_credit=arm_pair_credit, passed=all(t["pass_gate"] for t in tests.values()))


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    a=p.parse_args()
    dest=HERE/"decision.json"
    result=decision()
    if a.write:
        dest.write_text(json.dumps(result,indent=2)+"\n")
    else:
        assert json.loads(dest.read_text()) == result
    print(json.dumps(result, indent=2))
