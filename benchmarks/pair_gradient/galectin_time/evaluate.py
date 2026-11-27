"""Audit all raw score/cost records and rank experimentally separated compound pairs."""

import argparse
import csv
import hashlib
import itertools
import json
import re
import sys
import tarfile
from decimal import Decimal
from pathlib import Path
from statistics import median

from scipy.stats import kendalltau, spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from assay_transfer import RECEIVERS, SOURCE, save

FIELDS = (
    "hunt_evals",
    "second_evals",
    "explicit_evals",
    "hunt_calls",
    "second_calls",
    "explicit_calls",
)
EVALS = FIELDS[:3]


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def verify(work):
    budgets = {
        (r["receptor"], r["compound"]): int(r["calibrated_cap"])
        for r in read(work / "budgets.csv")
    }
    historic = {
        (r["receptor"], r["compound"], r["seed"], r["arm"]): r
        for r in read(HERE.parent / "galectins/cells.csv")
        if r["arm"] in ("parent", "linear")
    }
    cells, calibration = read(work / "validation.csv"), read(work / "calibration.csv")
    assert len(calibration) == 208 and len(cells) == 104 and len(budgets) == 52
    assert (
        json.loads((work / "run.json").read_text())["budget_sha256"]
        == hashlib.sha256((work / "budgets.csv").read_bytes()).hexdigest()
    )
    times = {
        (r["receptor"], r["compound"], r["arm"], r["cap"]): int(r["elapsed_ns"])
        for r in calibration
    }
    assert len(times) == 208
    nreplay = 0
    score_map = {}
    for r in calibration + cells:
        key, rec, c, arm, seed, cap = (
            r[f] for f in ("key", "receptor", "compound", "arm", "seed", "cap")
        )
        nominal = historic[rec, c, seed, arm]
        modes = re.findall(
            r"^REMARK VINA RESULT:\s+(\S+)",
            (work / "raw" / (key + ".pdbqt")).read_text(),
            re.M,
        )
        log = (work / "raw" / (key + ".log")).read_text()
        counts = {f: re.findall(r"\b" + f + r"=(\d+)", log) for f in FIELDS}
        if (
            len(modes) != int(r["modes"])
            or modes[0] != r["first_score"]
            or any(counts[f] != [r[f]] for f in FIELDS)
        ):
            raise ValueError("Raw outcome/counter disagreement: " + key)
        if cap == "12000":
            if Decimal(modes[0]) != Decimal(nominal["top_score"]) or any(
                r[f] != nominal[f] for f in FIELDS
            ):
                raise ValueError("Historical fixed-cap replay disagreement: " + key)
            nreplay += 1
        if seed == "98251":
            expected = 12000 if arm == "parent" else budgets[rec, c]
            if int(cap) != expected:
                raise ValueError("Validation cap changed")
            score_map[rec, int(c), arm] = Decimal(modes[0])
            if arm == "linear":
                score_map[rec, int(c), "linear_12k_history"] = Decimal(
                    nominal["top_score"]
                )
    assert len(score_map) == 156
    for (rec, c), got in budgets.items():
        t = times[rec, c, "parent", "12000"]
        lo, hi = (times[rec, c, "linear", x] for x in ("3000", "12000"))
        expected = (
            12000
            if t >= hi
            else 3000
            if t <= lo
            else int(3000 + Decimal(9000) * (t - lo) / (hi - lo))
        )
        if expected != got:
            raise ValueError("Calibration arithmetic mismatch")
    return cells, score_map, nreplay


def rank(work, scores):
    labels = {int(r["compound"]): r for r in read(SOURCE / "table1.csv")}
    with tarfile.open(SOURCE / "evidence.tar.xz") as archive:
        panel = {
            int(r["compound"]): r
            for r in csv.DictReader(
                archive.extractfile("panel.csv").read().decode().splitlines()
            )
        }
    summaries, pairs = [], []
    for target, recs in RECEIVERS.items():
        kd = {i: Decimal(r[target + "_Kd_uM"]) for i, r in labels.items()}
        assert kd == {i: Decimal(r[target + "_Kd_uM"]) for i, r in panel.items()}
        for state in (*recs, "consensus"):
            used = recs if state == "consensus" else (state,)
            for arm in ("parent", "linear", "linear_12k_history"):
                val = {
                    i: sum(scores[r, i, arm] for r in used) / len(used) for i in labels
                }
                for subset, ids in (
                    ("all", sorted(labels)),
                    ("triazoles", list(range(5, 17))),
                    ("equal_heavy", sorted(labels)),
                ):
                    edges = []
                    for a, b in itertools.combinations(sorted(ids, key=kd.get), 2):
                        if kd[b] < 3 * kd[a] or (
                            subset == "equal_heavy"
                            and panel[a]["heavy_atoms"] != panel[b]["heavy_atoms"]
                        ):
                            continue
                        gap = val[b] - val[a]
                        credit = float(gap > 0) if gap else 0.5
                        edges.append((a, b, credit))
                        pairs.append(
                            dict(
                                target=target,
                                state=state,
                                arm=arm,
                                subset=subset,
                                strong=a,
                                weak=b,
                                fold=float(kd[b] / kd[a]),
                                strong_score=str(val[a]),
                                weak_score=str(val[b]),
                                credit=credit,
                            )
                        )
                    x, y = [float(kd[i]) for i in ids], [float(val[i]) for i in ids]
                    summaries.append(
                        dict(
                            target=target,
                            state=state,
                            arm=arm,
                            subset=subset,
                            credit=sum(e[2] for e in edges),
                            pairs=len(edges),
                            spearman=""
                            if subset == "equal_heavy"
                            else float(spearmanr(x, y).statistic),
                            kendall_tau_b=""
                            if subset == "equal_heavy"
                            else float(kendalltau(x, y).statistic),
                        )
                    )
    save(work / "summary.csv", summaries)
    save(work / "pairs.csv", pairs)
    sensitivity = []
    for target in RECEIVERS:
        all_pairs = [
            r
            for r in pairs
            if r["target"] == target
            and r["state"] == "consensus"
            and r["subset"] == "all"
            and r["arm"] in ("linear", "parent")
        ]
        for omit in [None, *labels]:
            selected = [r for r in all_pairs if omit not in (r["strong"], r["weak"])]
            linear = [r for r in selected if r["arm"] == "linear"]
            parent = [r for r in selected if r["arm"] == "parent"]
            hc = [
                (
                    int(panel[r["strong"]]["heavy_atoms"]),
                    int(panel[r["weak"]]["heavy_atoms"]),
                )
                for r in parent
            ]
            sensitivity.append(
                {
                    "target": target,
                    "omit_compound": omit if omit else "none",
                    "pairs": len(parent),
                    "parent": sum(r["credit"] for r in parent),
                    "linear": sum(r["credit"] for r in linear),
                    "larger_ligand": sum(int(a > b) + 0.5 * int(a == b) for a, b in hc),
                }
            )
    save(work / "influence.csv", sensitivity)


def costs(work, cells):
    comp, groups = [], []
    index = {(r["receptor"], r["compound"], r["arm"]): r for r in cells}
    for target, recs in RECEIVERS.items():
        for rec, c in itertools.product(recs, map(str, range(5, 18))):
            a, b = (index[rec, c, arm] for arm in ("parent", "linear"))
            tp, tl = (int(r["elapsed_ns"]) / 1e9 for r in (a, b))
            ep, el = (sum(int(r[f]) for f in EVALS) for r in (a, b))
            comp.append(
                dict(
                    target=target,
                    receptor=rec,
                    compound=c,
                    linear_cap=int(b["cap"]),
                    parent_s=tp,
                    linear_s=tl,
                    ratio=tl / tp,
                    parent_evals=ep,
                    linear_evals=el,
                )
            )
    for name in [*RECEIVERS, *sum(RECEIVERS.values(), ()), "all"]:
        rows = [r for r in comp if name in (r["target"], r["receptor"], "all")]
        sums = {
            k: sum(r[k] for r in rows)
            for k in ("parent_s", "linear_s", "parent_evals", "linear_evals")
        }
        ratios = [r["ratio"] for r in rows]
        groups.append(
            dict(
                group=name,
                cells=len(rows),
                **sums,
                total_ratio=sums["linear_s"] / sums["parent_s"],
                eval_ratio=sums["linear_evals"] / sums["parent_evals"],
                median_cell_ratio=median(ratios),
                min_cell_ratio=min(ratios),
                max_cell_ratio=max(ratios),
                within_5pp=sum(0.95 <= v <= 1.05 for v in ratios),
                no_slower=sum(v <= 1 for v in ratios),
            )
        )
    save(work / "paired_costs.csv", comp)
    save(work / "costs.csv", groups)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    work = parser.parse_args().work
    cells, scores, replays = verify(work)
    rank(work, scores)
    costs(work, cells)
    print(
        f"Audited 312 raw pose/CLI outputs; {replays} fixed-cap scores and six counters match the historical run."
    )
