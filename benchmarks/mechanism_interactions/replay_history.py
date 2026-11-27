"""Recompute selected historical interactions from per-observation records (no docking)."""

import argparse
import csv
import itertools
import lzma
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from scipy.stats import kendalltau, spearmanr


HERE = Path(__file__).resolve().parent
BENCH = HERE.parent


def records(path):
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "rt", newline="") as stream:
        yield from csv.DictReader(stream)


def save(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: format(v, ".10g") if isinstance(v, float) else v for k, v in row.items()})


def population_ratio(ph1, ph2, ligand_pka, *, ligand=True, receptor=True):
    # Ratio of free binding-competent fractions at reference pH to test pH.
    lr = (1 + 10 ** (ligand_pka - ph2)) / (1 + 10 ** (ligand_pka - ph1))
    rr = (1 + 10 ** (ph2 - 7.01)) / (1 + 10 ** (ph1 - 7.01))
    return (lr if ligand else 1.0) * (rr if receptor else 1.0)


def carbonic_anhydrase(out):
    src = BENCH / "proton_linkage"
    series = list(records(src / "caix.csv"))
    predictions, performance, counterexamples = [], [], []
    for arm, ligand, receptor in (("reference", False, False), ("ligand", True, False),
                                  ("receptor", False, True), ("both", True, True)):
        errors = []
        for row in series:
            kd = float(row["kd_ph5_nm"])
            pred = float(row["kd_ph7_nm"]) * population_ratio(7, 5, float(row["ligand_pka"]),
                                                               ligand=ligand, receptor=receptor)
            err = math.log10(pred / kd)
            errors.append(abs(err))
            predictions.append(dict(ligand=row["ligand"], arm=arm, kd_pH5_nM=kd, predicted_nM=pred,
                                    log10_error=err))
        performance.append(dict(arm=arm, ligands=len(errors), mae_log10=mean(errors)))
    pka13 = next(float(row["ligand_pka"]) for row in series if row["ligand"] == "13")
    for row in records(src / "caix_itc.csv"):
        ref, test = float(row["reference_ph"]), float(row["test_ph"])
        measured, initial = float(row["test_kd_nm"]), float(row["reference_kd_nm"])
        counterexamples.append(dict(reference=row["reference"], reference_pH=ref, test_pH=test,
                                    measured_nM=measured, predicted_nM=initial*population_ratio(ref,test,pka13),
                                    observed_ratio=measured/initial, max_ratio=10**abs(test-ref)))
    save(out / "caix_predictions.csv", predictions)
    save(out / "caix_mae.csv", performance)
    save(out / "caix_alkaline.csv", counterexamples)


def rigid_hinge(out):
    src = BENCH / "affinity_reference_ensemble" / "results"
    paths = [src / "primary_historical_cells.csv.xz", src / "s6k_cells.csv.xz"]
    scores = defaultdict(lambda: defaultdict(list))
    endpoints = {}
    observed = set()
    for row in itertools.chain.from_iterable(map(records, paths)):
        if int(row["nodes"]) != 1024:
            continue
        target, ligand, receptor = row["target"], row["ligand"], row["state"]
        key = (target, ligand, receptor, row["seed"])
        if key in observed:
            raise ValueError(f"Repeated ligand/receptor/seed: {key}")
        observed.add(key)
        kd = float(row["kd_nM"])
        if endpoints.setdefault((target, ligand), kd) != kd:
            raise ValueError(f"Inconsistent Kd for {target} {ligand}")
        native, rigid, hinge, joint = map(float, (row[k] for k in ("native", "rigid", "hinge", "joint")))
        arms = dict(native=native, rigid=native+rigid, hinge=native+hinge,
                    separate=native+rigid+hinge, joint=native+joint, size=-int(row["n_heavy"]))
        for arm, score in arms.items():
            scores[(target, arm, ligand)][receptor].append(score)

    ligand_scores = defaultdict(dict)
    for (target, arm, ligand), receivers in scores.items():
        ligand_scores[(target, arm)][ligand] = mean(median(values) for values in receivers.values())

    rows = []
    for (target, arm), by_ligand in sorted(ligand_scores.items()):
        order = sorted(by_ligand, key=lambda l: endpoints[(target, l)])
        kd = [endpoints[(target, l)] for l in order]
        score = [by_ligand[l] for l in order]
        credits, errors = [], []
        for i, j in itertools.combinations(range(len(order)), 2):
            if kd[j] < 3 * kd[i]:
                continue
            gap = score[j] - score[i]
            credits.append(1 if gap > 0 else (0.5 if gap == 0 else 0))
            errors.append(abs(gap - 0.00198720425864083 * 298.15 * math.log(kd[j] / kd[i])))
        rows.append(dict(target=target, arm=arm, ligands=len(order), pairs_ge3=len(credits),
                         pair_credit=sum(credits), spearman=spearmanr(score, kd).statistic,
                         kendall_b=kendalltau(score, kd).statistic,
                         gap_mae_kcal_mol=mean(errors) if arm != "size" else ""))
    save(out / "rigid_hinge.csv", rows)


def archive(out):
    src = BENCH / "archive_overlap" / "transfer"
    by_arm = defaultdict(dict)
    for directory in (src, src / "mechanism" / "factorial"):
        modes = defaultdict(dict)
        for row in records(directory / "modes.csv"):
            modes[row["key"]][int(row["mode"])] = float(row["rmsd_A"])
        for row in records(directory / "cells.csv"):
            arm = {"candidate": "both", "parent": "baseline", "save": "save", "merge": "merge"}.get(row["arm"])
            if arm is None:
                continue
            rs = modes[row["key"]]
            if sorted(rs) != list(range(1, int(row["modes"]) + 1)):
                raise ValueError(f"Incomplete ranked pool for {row['key']}")
            case = (row["target"], row["query"], row["receiver"], row["seed"])
            by_arm[arm][case] = (rs[1], min(r for k, r in rs.items() if k <= 5), min(rs.values()),
                                 sum(int(row[k]) for k in ("hunt_evals", "second_evals", "explicit_evals")))
    cases = set(by_arm["baseline"])
    if any(set(pool) != cases for pool in by_arm.values()):
        raise ValueError("Archive arms lack matched structural conditions")
    rows = []
    for arm, pool in sorted(by_arm.items()):
        for target in ["all", *sorted({c[0] for c in cases})]:
            treatment = [v for c, v in pool.items() if target in ("all", c[0])]
            top, best, oracle, evals = map(list, zip(*treatment))
            rows.append(dict(target=target, arm=arm, cells=len(top), mean_top1_A=mean(top),
                             median_top1_A=median(top), mean_best5_A=mean(best), mean_returned_oracle_A=mean(oracle),
                             top1_le2=sum(r <= 2 for r in top), best5_le2=sum(r <= 2 for r in best),
                             oracle_le2=sum(r <= 2 for r in oracle), evaluations=sum(evals)))
    save(out / "archive.csv", rows)


def gba1(out):
    src = BENCH / "ch_pi_transfer" / "gba1" / "results"
    by_budget = {}
    for cap, directory in ((12000, src), (24000, src / "cap24000")):
        pools = defaultdict(list)
        for r in records(directory / "modes.csv"):
            if (r["state"], r["engine"]) == ("protonated", "pro"):
                pools[(r["query"], r["receptor"], r["seed"])].append(r)
        by_budget[cap] = {}
        for key, rs in pools.items():
            ranks = sorted(int(r["rank"]) for r in rs)
            if ranks != list(range(1,10)):
                raise ValueError(f"Expected nine unique modes for {key} {cap}: {ranks}")
            top = min(rs, key=lambda r: int(r["rank"]))
            chpi = min(rs, key=lambda r: (float(r["chpi"]), int(r["rank"])))
            by_budget[cap][key] = (float(top["rmsd"]), float(chpi["rmsd"]),
                                    min(float(r["rmsd"]) for r in rs))
    if by_budget[12000].keys() != by_budget[24000].keys():
        raise ValueError("GBA1 search budgets lack matched conditions")
    rows = []
    for cap, vals in by_budget.items():
        native, chpi, oracle = map(list, zip(*vals.values()))
        delta = [b-a for a,b in zip(native, chpi)]
        rows.append(dict(cap=cap, cells=len(delta), mean_native_A=mean(native), mean_chpi_A=mean(chpi),
                         median_native_A=median(native), median_chpi_A=median(chpi), mean_pool_best_A=mean(oracle),
                         mean_delta_A=mean(delta), median_delta_A=median(delta), better=sum(x<0 for x in delta),
                         worse=sum(x>0 for x in delta), equal=sum(x==0 for x in delta)))
    save(out / "gba1.csv", rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "results")
    out = parser.parse_args().output
    out.mkdir(parents=True, exist_ok=True)
    for replay in (carbonic_anhydrase, rigid_hinge, archive, gba1):
        replay(out)
