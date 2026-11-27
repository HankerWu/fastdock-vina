"""Check raw scores and work counters, then calculate paired within-assay metrics."""

import argparse
import csv
import itertools
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from scipy.stats import kendalltau, spearmanr

SEED = 458921


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w") as stream:
        output = csv.DictWriter(stream, list(rows[0]), lineterminator="\n")
        output.writeheader()
        output.writerows(rows)


def credit(gap):
    return float(gap > 0) + 0.5 * (gap == 0)


def analyze(work, output, cost_check=False):
    study = "cost_check" if cost_check else "experiments"
    seed = 703109 if cost_check else SEED
    output.mkdir(parents=True, exist_ok=True)
    runs = read(work / study / "runs.csv")
    panel = {(r["target"], r["cid"]): r for r in read(work / "panel.csv")}
    groups = defaultdict(dict)
    for r in runs:
        stem = work / study / "raw" / r["key"]
        ranks = re.findall(
            r"^REMARK VINA RESULT:\s+([-.\d]+)",
            stem.with_suffix(".pdbqt").read_text(),
            re.MULTILINE,
        )
        assert float(ranks[0]) == float(r["first_score"]) and len(ranks) == int(
            r["modes"]
        )
        log = stem.with_suffix(
            ".counts" if r["arm"] == "upstream" else ".log"
        ).read_text()
        for k in (
            "hunt_calls",
            "second_calls",
            "explicit_calls",
            "hunt_evals",
            "second_evals",
            "explicit_evals",
        ):
            assert int(re.search(r"\b" + k + r"=(\d+)", log)[1]) == int(r[k])
        assert sum(
            int(r[k]) for k in ("hunt_evals", "second_evals", "explicit_evals")
        ) == int(r["joint"])
        r["score"] = float(ranks[0])
        r["kd"] = float(panel[r["panel"], r["compound"]]["kd_nm"])
        key = r["panel"], int(r["cap"]), r["arm"]
        assert r["compound"] not in groups[key] and int(r["seed"]) == seed
        groups[key][r["compound"]] = r
    expected = set(
        itertools.product(
            {t for t, _ in panel},
            (10800, 12000) if cost_check else (3000, 12000),
            ("parent", "halving", "quadratic")
            if cost_check
            else ("upstream", "parent", "halving", "quadratic"),
        )
    )
    assert set(groups) == expected
    summary, edges, drops, work_rows = [], [], [], []
    for (target, cap, arm), cells in sorted(groups.items()):
        base = {"target": target, "cap": cap, "arm": arm}
        assert cells.keys() == {cid for t, cid in panel if t == target}
        ligs = sorted(cells, key=lambda cid: cells[cid]["kd"])
        subset = []
        for a, b in itertools.combinations(ligs, 2):
            strong, weak = cells[a], cells[b]
            ratio = weak["kd"] / strong["kd"]
            if ratio < 3:
                continue
            sa, wb = panel[target, a], panel[target, b]
            gap = weak["score"] - strong["score"]
            ref = 0.00198720425864 * 298.15 * math.log(ratio)
            row = dict(
                **base,
                strong=a,
                weak=b,
                kd_ratio=ratio,
                gap=round(gap, 6),
                ref_gap_298K=ref,
                abs_error=abs(gap - ref),
                credit=credit(gap),
                size_credit=credit(int(sa["heavy"]) - int(wb["heavy"])),
                same_size_charge=int(
                    (sa["heavy"], sa["charge"]) == (wb["heavy"], wb["charge"])
                ),
                same_scaffold=int(sa["scaffold"] == wb["scaffold"]),
                original=int(sa["old_panel"] == wb["old_panel"] == "1"),
                without_cocrystal=int(sa["cocrystal"] == wb["cocrystal"] == "0"),
            )
            subset.append(row)
        edges += subset
        for selected in (
            "all",
            "without_cocrystal",
            "same_size_charge",
            "same_scaffold",
            "cross_scaffold",
            "original",
        ):
            ss = [
                p
                for p in subset
                if selected == "all"
                or (
                    not p["same_scaffold"]
                    if selected == "cross_scaffold"
                    else p[selected]
                )
            ]
            summary.append(
                dict(
                    **base,
                    subset=selected,
                    pairs=len(ss),
                    correct=sum(p["credit"] for p in ss),
                    size=sum(p["size_credit"] for p in ss),
                    gap_mae=mean(p["abs_error"] for p in ss) if ss else "",
                )
            )
        for cid in ligs:
            leave = [p for p in subset if cid not in (p["strong"], p["weak"])]
            drops.append(
                dict(
                    **base,
                    removed=cid,
                    n=len(leave),
                    correct=sum(p["credit"] for p in leave),
                )
            )
        kd = [cells[c]["kd"] for c in ligs]
        predicted = [cells[c]["score"] for c in ligs]
        work_rows.append(
            dict(
                **base,
                n=len(ligs),
                spearman=spearmanr(kd, predicted).statistic,
                kendall_b=kendalltau(kd, predicted).statistic,
                total_s=sum(int(r["elapsed_ns"]) for r in cells.values()) / 1e9,
                median_ms=median(int(r["elapsed_ns"]) for r in cells.values()) / 1e6,
                joint=sum(int(r["joint"]) for r in cells.values()),
                local_calls=sum(
                    sum(
                        int(r[k])
                        for k in ("hunt_calls", "second_calls", "explicit_calls")
                    )
                    for r in cells.values()
                ),
            )
        )
    for name, rows in (
        ("summary", summary),
        ("pairs", edges),
        ("leave_one_out", drops),
        ("cost_ranks", work_rows),
    ):
        write(output / (name + ".csv"), rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--cost-check", action="store_true")
    args = parser.parse_args()
    analyze(args.work, args.out, args.cost_check)
