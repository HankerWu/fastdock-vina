"""Re-score three measured assays without changing any coordinates or unbound reference."""
import argparse
from collections import defaultdict
import csv
from decimal import Decimal
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from statistics import mean
import sys
import time

from scipy.stats import kendalltau, spearmanr

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "hydroxyl_partition"))
from scorer import Receptor, score_pose
sys.path.insert(0, str(BENCH / "terminal_rotor_score"))
from assay import WEIGHT, aggregate, assay_inputs, credit, read, save

ARMS = ("native", "fixed", "minimum", "uniform", "ensemble", "rotor_only")


def score():
    records = {}
    for r in read(BENCH / "terminal_rotor_score/results/poses.csv"):
        records[r["target"], r["ligand"], r["receptor"], int(r["seed"])] = r
    rows, costs = [], []

    @lru_cache(maxsize=None)
    def prepare(content):
        start = time.process_time()
        result = Receptor(content.decode())
        costs.append({"stage": "receptor", "cpu_seconds": time.process_time() - start})
        return result

    for target, ligand, measured, receptor, seed, protein, topology, modes in assay_inputs():
        rec = prepare(protein)
        source = records.pop((target, ligand, receptor, seed))
        native = float(source["fixed"])
        scale = 1 / (1 + WEIGHT * float(source["fixed_rotors"]))
        first = modes.decode().split("ENDMDL", 1)[0]
        if abs(float(source["saved_inter"]) * scale - native) > 1.e-10:
            raise ValueError("Inconsistent native scoring reference")
        if float(first.split("REMARK VINA RESULT:")[1].split()[0]) != float(source["saved_score"]):
            raise ValueError("Incorrect archived model")
        start = time.process_time()
        terms = score_pose(rec, topology.decode(), first)
        costs.append({"stage": "first_pose", "cpu_seconds": time.process_time() - start})
        rows.append(dict(target=target, ligand=ligand, measurement=measured, receptor=receptor,
                         seed=seed, scale=scale, native=native, **terms,
                         **{arm: native + scale * terms["raw_" + arm] for arm in ARMS[1:]}))
    if records or len(rows) != 234:
        raise ValueError("Incomplete 234-pose panel")
    return rows, costs


def summarize(rows):
    summaries, edges = [], []
    for target in sorted({r["target"] for r in rows}):
        panel = [r for r in rows if r["target"] == target]
        for rec in ["mean", *sorted({r["receptor"] for r in panel})]:
            cells = [r for r in panel if rec == "mean" or r["receptor"] == rec]
            exp = {r["ligand"]: Decimal(str(r["measurement"])) for r in cells}
            keys = sorted(exp, key=lambda key: (exp[key], key))
            scores = {arm: aggregate(cells, arm) for arm in ARMS}
            pairs = []
            for strong, weak in combinations(keys, 2):
                if exp[weak] < 3 * exp[strong]:
                    continue
                row = dict(target=target, receptor=rec, strong=strong, weak=weak,
                           measured_fold=float(exp[weak] / exp[strong]))
                for arm in ARMS:
                    row[arm] = scores[arm][weak] - scores[arm][strong]
                edges.append(row)
                pairs.append(row)
            for arm in ARMS:
                delta = [credit(p[arm]) - credit(p["native"]) for p in pairs]
                leave = [sum(d for p, d in zip(pairs, delta)
                             if key not in (p["strong"], p["weak"])) for key in keys]
                obs = [float(exp[k]) for k in keys]
                pred = [scores[arm][k] for k in keys]
                summaries.append(dict(target=target, receptor=rec, arm=arm, ligands=len(keys),
                    cells=len(cells), mobile_cells=sum(bool(r["rotors"]) for r in cells),
                    edges=len(pairs), correct=sum(credit(p[arm]) for p in pairs),
                    gains=sum(max(d, 0.) for d in delta), losses=-sum(min(d, 0.) for d in delta),
                    spearman=float(spearmanr(obs, pred).statistic), kendall=float(kendalltau(obs, pred).statistic),
                    min_leave_one_ligand_delta=min(leave), max_leave_one_ligand_delta=max(leave),
                    mean_abs_change=mean(abs(scores[arm][k] - scores["native"][k]) for k in keys)))
    return summaries, edges


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rows, costs = score()
    sums, edges = summarize(rows)
    for name, data in (("affinity_poses", rows), ("affinity_summary", sums), ("affinity_cost", costs)):
        save(args.output / (name + ".csv"), data)
    import lzma
    with lzma.open(args.output / "affinity_pairs.csv.xz", "wt", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=edges[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(edges)


if __name__ == "__main__":
    main()
