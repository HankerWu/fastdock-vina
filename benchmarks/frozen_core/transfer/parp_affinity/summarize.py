"""Check native C++ rotor annotation vs fixed-geometry Python and compare within-assay scores."""

import argparse
import hashlib
import itertools
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from statistics import median

from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "affinity"))
import test_entropy as check

SEEDS = (98207, 98213)
STATES = ("4UND", "7KK4")


def extraction():
    r = ET.parse(HERE / "table1.xml")
    labels = {}
    for row in r.findall(".//tbody/tr"):
        fields = ["".join(c.itertext()).strip() for c in row]
        value, error = map(float, re.split(r"\s*±\s*", fields[4]))
        labels[fields[0].lower()] = (value, error)
    assert len(labels) == 4
    return labels


def analyze(work, parent, out):
    labels = extraction()
    cc = {r["compound"]: r for r in check.read_csv(HERE / "native_rotors.csv")}
    mmff = {
        (v["compound"], v["receptor"], int(v["seed"])): v
        for v in check.read_csv(HERE / "mmff.csv")
    }
    receptors = {
        s: cKDTree(
            [
                a[2]
                for a in map(
                    check.atom,
                    check.atom_lines((parent / s / "receptor.pdbqt").read_text()),
                )
                if check.heavy(a)
            ]
        )
        for s in STATES
    }
    selected = []
    modes = []
    for l, s, k in itertools.product(labels, STATES, SEEDS):
        pool = (work / "pools" / f"{l}-{s}-{k}.pdbqt").read_text()
        block = pool.split("ENDMDL")[0]
        topo = check.ligand_graph((work / (l + ".pdbqt")).read_text())
        for i, mode in enumerate((m for m in pool.split("ENDMDL") if m.strip()), 1):
            values = check.pose(mode, receptors[s], topo, float(cc[l]["num_tors"]))
            (remark,) = [
                line
                for line in mode.splitlines()
                if line.startswith("REMARK VINA CONTACT ROTOR SCORE:")
            ]
            energy, n, kinds = remark.split(":")[1].split()
            diff = abs(float(energy) - values["contact"])
            assert (int(n), int(kinds)) == (
                values["constrained"],
                values["eligible_rotors"],
            ) and diff < 0.00101
            modes.append(
                dict(
                    compound=l,
                    receptor=s,
                    seed=k,
                    model=i,
                    cpp=float(energy),
                    delta=diff,
                    **values,
                )
            )
        terms = check.pose(block, receptors[s], topo, float(cc[l]["num_tors"]))
        (line,) = [
            z
            for z in block.splitlines()
            if z.startswith("REMARK VINA CONTACT ROTOR SCORE:")
        ]
        native_contact, constrained, eligible = line.split(":")[1].split()
        assert (
            int(constrained) == terms["constrained"]
            and int(eligible) == terms["eligible_rotors"]
        )
        assert abs(float(native_contact) - terms["contact"]) < 0.00101
        mm = mmff[l, s, k]
        assert abs(float(mm["native_score"]) - terms["native"]) < 0.00101
        selected.append(
            dict(
                compound=l,
                receptor=s,
                seed=k,
                kd_nM=labels[l][0],
                reported_plusminus_nM=labels[l][1],
                pool_sha256=hashlib.sha256(pool.encode()).hexdigest(),
                mmff=float(mm["unit_score"]),
                native=terms["native"],
                raw=terms["raw"],
                uniform=terms["uniform"],
                contact=float(native_contact),
                constrained=terms["constrained"],
                eligible=terms["eligible_rotors"],
                num_tors=terms["num_tors"],
            )
        )
    check.write(out / "rotor_replay.csv", modes)
    measured = {r["compound"]: r["kd_nM"] for r in selected}
    preds, pairs, summary = [], [], []
    for policy in [*STATES, "mean_receptors"]:
        states = STATES if policy == "mean_receptors" else [policy]
        p = {}
        for l in measured:
            p[l] = {}
            for method in ("native", "raw", "uniform", "contact", "mmff"):
                p[l][method] = sum(
                    median(
                        r[method]
                        for r in selected
                        if r["compound"] == l and r["receptor"] == s
                    )
                    for s in states
                ) / len(states)
            preds.append(dict(receptor=policy, compound=l, kd_nM=measured[l], **p[l]))
        for method in ("native", "raw", "uniform", "contact", "mmff"):
            group = []
            for a, b in itertools.combinations(sorted(measured, key=measured.get), 2):
                ratio = measured[b] / measured[a]
                delta = round(p[b][method] - p[a][method], 10)
                row = dict(
                    receptor=policy,
                    method=method,
                    stronger=a,
                    weaker=b,
                    fold=ratio,
                    eligible=int(ratio >= 3),
                    delta=delta,
                    credit=1 if delta > 0 else 0 if delta < 0 else 0.5,
                )
                pairs.append(row)
                if ratio >= 3:
                    group.append(row)
            assert len(group) == 4
            summary.append(
                dict(
                    receptor=policy,
                    method=method,
                    pairs=len(group),
                    correct=sum(row["credit"] for row in group),
                )
            )
    for name, rows in [
        ("selected", selected),
        ("predictions", preds),
        ("pairs", pairs),
        ("summary", summary),
    ]:
        check.write(out / (name + ".csv"), rows)
    print(summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("parent", type=Path)
    p.add_argument("--out", type=Path, default=HERE)
    a = p.parse_args()
    analyze(a.work, a.parent, a.out)
