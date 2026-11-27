"""Joint apo-population / common ligand-reference impossibility certificates.

Run from this directory or any working directory. All comparisons use the
frozen native scores, not fitted binding free energies. --frozen is a replay
of the small archived state-median table when the original archives are absent.
"""
import argparse
import csv
import itertools
import json
import math
import sys
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
from statistics import median

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
PANELS = [("Gal1", "Gal3"), ("AURKA", "CHEK1", "S6K1")]
STATES = {"Gal1": ("1GZW", "4Q26"), "Gal3": ("1KJL", "3T1L"),
          "AURKA": ("2X6E", "2J4Z"), "CHEK1": ("2E9P", "2YWP"), "S6K1": ("3A60",)}
RT = 0.00198720425864083 * 298.15


def table(path):
    with path.open(newline="") as stream:
        yield from csv.DictReader(stream)


def load_states():
    # Only exact same-assay Kd entries retained in the historical 500-cell panel.
    primary = list(table(BENCH / "affinity_reference_ensemble/results/primary_ligands.csv"))
    allowed = {(r["target"], int(r["ligand"])) for r in primary}
    samples = defaultdict(dict)
    for r in table(BENCH / "galectin_affinity/results/readback.csv"):
        if r["arm"] == "native":
            t = next(t for t in ("Gal1", "Gal3") if r["state"] in STATES[t])
            k = (t, int(r["compound"]), r["state"])
            if r["seed"] in samples[k]:
                raise ValueError("duplicate Galectin seed")
            samples[k][r["seed"]] = F(r["saved"])
    for r in table(BENCH / "torsion_kinase_affinity/cells.csv"):
        t, i = r["target"], int(r["cid"])
        if r["arm"] == "native" and (t, i) in allowed:
            k = (t, i, r["state"])
            if r["seed"] in samples[k]:
                raise ValueError("duplicate kinase seed")
            samples[k][r["seed"]] = F(r["pool_score"])
    data = defaultdict(dict)
    for r in primary:
        t, i = r["target"], int(r["ligand"])
        if t == "S6K1":  # Only one archived receptor; native is its seed median.
            s = (F(r["native"]),)
        else:
            ss = [samples[t, i, state] for state in STATES[t]]
            if any(len(v) != 4 for v in ss):
                raise ValueError(f"incomplete four-seed state: {t}, {i}")
            s = tuple(median(v.values()) for v in ss)
        if abs(sum(s) / len(s) - F(r["native"])) > F("0.00000001"):
            raise ValueError(f"published native median mismatch: {t}, {i}")
        if i in data[t]:
            raise ValueError("duplicate target-ligand entry")
        data[t][i] = (F(r["kd_nM"]), s)
    return data


def load_frozen(path):
    data = defaultdict(dict)
    for r in table(path):
        t, i = r["target"], int(r["ligand"])
        s = tuple(F(r[k]) for k in ("score_1", "score_2") if r[k])
        if len(s) != len(STATES[t]) or i in data[t]:
            raise ValueError("invalid state table")
        data[t][i] = (F(r["kd_nM"]), s)
    return data


def gaps(data, t, i, j):
    return tuple(a - b for a, b in zip(data[t][i][1], data[t][j][1]))


def analyze(data, output):
    # Reuse exact maximum-cycle solver and independent HiGHS LP certificate.
    sys.path.insert(0, str(HERE.parent))
    from common_reference.cycles import maximum_cycle_mean, certify

    if {t: len(v) for t, v in data.items()} != {"Gal1": 13, "Gal3": 13, "AURKA": 19, "CHEK1": 14, "S6K1": 14}:
        raise ValueError("unexpected panel coverage")
    result = {"states": STATES, "pairs": [], "impossible_reversals": [], "graphs": []}
    for targets in PANELS:
        removed = set()
        for t, u in itertools.combinations(targets, 2):
            n = rev = bad = strong_bad = 0
            errors = []
            for i, j in itertools.combinations(sorted(data[t].keys() & data[u].keys()), 2):
                x, y = data[t][i][0], data[t][j][0]
                z, w = data[u][i][0], data[u][j][0]
                a, b = gaps(data, t, i, j), gaps(data, u, i, j)
                lo, hi = min(a)-max(b), max(a)-min(b)
                ratio = x*w/(y*z)
                obs = RT*math.log(float(ratio))
                wrong = (ratio > 1 and hi <= 0) or (ratio < 1 and lo >= 0)
                errors.append(max(float(lo)-obs, obs-float(hi), 0))
                n += 1
                if max(ratio, 1/ratio) >= 3:
                    strong_bad += wrong
                if (x >= 3*y and w >= 3*z) or (y >= 3*x and z >= 3*w):
                    rev += 1
                    if wrong:
                        bad += 1
                        removed.add((i, j))
                        result["impossible_reversals"].append(dict(targets=[t,u], ligands=[i,j], kd_nM=[float(v) for v in (x,y,z,w)], experimental=obs, lower=float(lo), upper=float(hi), exact_lower=str(lo), exact_upper=str(hi)))
            result["pairs"].append(dict(targets=[t,u], pairs=n, reversals=rev, impossible=bad, large_exchange_wrong=strong_bad, interval_mae_floor=sum(errors)/n))
        edges = [(i,j,t,min(gaps(data,t,i,j))) for t in targets for i,j in itertools.permutations(sorted(data[t]),2) if data[t][j][0] >= 3*data[t][i][0]]
        for name, es in [("all", edges), ("without_impossible_reversals", [e for e in edges if tuple(sorted(e[:2])) not in removed])]:
            mean, cycle = maximum_cycle_mean(es)
            lp, _ = certify(es, mean)
            if mean is not None and abs(lp-float(mean)) > 1e-8:
                raise ArithmeticError("LP disagrees with exact cycle")
            if cycle and sum(e[3] for e in cycle)/len(cycle) != mean:
                raise ArithmeticError("invalid witness")
            result["graphs"].append(dict(targets=targets, graph=name, edges=len(es), lower_bound=float(mean) if mean is not None else None, exact_lower_bound=str(mean) if mean is not None else None, universal_conflict=mean is not None and mean >= 0, witness=[dict(stronger=i,weaker=j,target=t,lower=float(w),exact_lower=str(w)) for i,j,t,w in cycle]))
    output.write_text(json.dumps(result, indent=2)+"\n")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frozen", type=Path, help="replay state medians instead of upstream raw CSVs")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    data = load_frozen(a.frozen) if a.frozen else load_states()
    analyze(data, a.output)
