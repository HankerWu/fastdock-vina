"""Reconcile experimental labels, top-ranked energies and paired search budgets."""
import argparse
import csv
import itertools
import re
from collections import defaultdict
from pathlib import Path
from statistics import median

from scipy.stats import spearmanr

from run import ARMS, SEEDS, STATES


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def summarize(work, seeds):
    with (work / "panel.csv").open() as stream:
        panel = list(csv.DictReader(stream))
    cells, predictions, pairs, summaries = [], [], [], []
    for row in panel:
        for state in STATES[row["target"]]:
            scores = defaultdict(list)
            for seed in seeds:
                provenance = []
                for arm in ARMS:
                    stem = work / "pools" / f"{state}_{row['compound']}_{seed}_{arm}"
                    logs = stem.with_suffix(".log").read_text().splitlines()
                    pdbqt = stem.with_suffix(".pdbqt").read_text().splitlines()
                    assert pdbqt[0] == "MODEL 1" and pdbqt[-1] == "ENDMDL"
                    score = float(pdbqt[1].split()[3])
                    matches = [line for line in logs if line.startswith(("Search trajectory:", "Hunt chain:", "Pocket start:"))]
                    assert len(matches) == 6
                    counters, = [line for line in logs if line.startswith("Search counters:")]
                    explicit = dict(re.findall(r"(explicit_\w+|late_prescore)=(\d+)", counters))
                    provenance.append((matches, explicit["explicit_inputs"], explicit["late_prescore"]))
                    scores[arm].append(score)
                    cells.append(dict(target=row["target"], assay=row["assay"], receptor=state, compound=row["compound"],
                                      seed=seed, arm=arm, top1_score=score, modes=sum(l.startswith("MODEL ") for l in pdbqt),
                                      **explicit))
                if provenance[0] != provenance[1]:
                    raise ValueError(f"Search differed across arms: {state}, {row['compound']}, {seed}")
            predictions.append(dict(target=row["target"], assay=row["assay"], receptor=state,
                                    compound=row["compound"], endpoint=row["endpoint"], value_nM=float(row["value_nM"]),
                                    **{arm: median(v) for arm, v in scores.items()}))
    grouped = defaultdict(list)
    for row in predictions:
        grouped[row["target"], row["assay"], row["receptor"]].append(row)
    for (target, assay, receptor), group in grouped.items():
        subset = []
        for left, right in itertools.combinations(group, 2):
            strong, weak = sorted((left, right), key=lambda r: r["value_nM"])
            ratio = weak["value_nM"] / strong["value_nM"]
            if ratio == 1:
                continue
            case = dict(target=target, assay=assay, receptor=receptor, stronger=strong["compound"],
                        weaker=weak["compound"], fold=ratio)
            for arm in ARMS:
                delta = weak[arm] - strong[arm]
                case[arm + "_delta"] = round(delta, 6)
                case[arm + "_correct"] = 1.0 if delta > 0 else 0.5 if delta == 0 else 0.0
            subset.append(case)
        pairs.extend(subset)
        for threshold in (1, 3, 10):
            selected = [pair for pair in subset if pair["fold"] >= threshold]
            for arm in ARMS:
                summaries.append(dict(target=target, assay=assay, receptor=receptor, arm=arm,
                                      compounds=len(group), fold_min=threshold, pairs=len(selected),
                                      correct=sum(r[arm + "_correct"] for r in selected),
                                      spearman=round(spearmanr([r["value_nM"] for r in group],
                                                               [r[arm] for r in group]).statistic, 5)))
    for filename, rows in (("cells", cells), ("predictions", predictions), ("pairs", pairs), ("summary", summaries)):
        write(work / (filename + ".csv"), rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--seeds", type=int, nargs="+", default=SEEDS)
    args = p.parse_args()
    summarize(args.work, args.seeds)
