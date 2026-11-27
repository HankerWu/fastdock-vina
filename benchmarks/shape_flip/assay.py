"""Redock all 47 members of the archived Gal3 Kd assay. Report receptor-stratified results."""

import argparse
import hashlib
import itertools
import os
import random
import re
import subprocess
import sys
import tempfile
from decimal import Decimal
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau, spearmanr

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH / "accepted_retention" / "affinity"))
from run import cells, dock, read, write


def main(a):
    os.sched_setaffinity(0, {a.cpu})
    path = a.work.resolve()
    (path / "raw").mkdir(parents=True, exist_ok=True)
    entries = cells(path, "Gal3")
    random.Random(617).shuffle(entries)
    arms = ["parent", "always", "adaptive", "upstream"]
    out = []
    labels = {
        r["compound"]: r
        for r in read(BENCH / "accepted_retention/affinity/inputs.csv")
        if r["panel"] == "Gal3"
    }
    for j, cell in enumerate([] if a.replay_scores else entries):
        for arm in arms[j % 4 :] + arms[: j % 4]:
            out.append(
                dock(cell, path, getattr(a, arm), arm, 391781, 3000, arm != "upstream")
            )
    if a.replay_scores:
        out = read(path / "runs.csv")
    else:
        write(path / "runs.csv", out)
    counted = []
    for cell in entries:
        row = dock(
            cell, path, a.upstream_counts, "upstream", 391781, 3000, False, diag=True
        )
        key = row["key"]
        replay = path / "raw" / (key + ".replay")
        if (
            replay.read_bytes() != (path / "raw" / (key + ".pdbqt")).read_bytes()
            or not row["joint"]
        ):
            raise ValueError(key + " counting replay mismatch")
        replay.unlink()
        counted.append(row)
    write(path / "upstream_counts.csv", counted)
    summary, pair_details = [], []
    for arm, rec in itertools.product(arms, ("1KJL", "3T1L", "mean")):
        scores = {}
        for c in labels:
            rows = [
                r
                for r in out
                if r["arm"] == arm
                and r["compound"] == c
                and (rec == "mean" or r["receptor"] == rec)
            ]
            scores[c] = sum(Decimal(r["first_score"]) for r in rows) / len(rows)
        kd = {c: Decimal(r["k"]) for c, r in labels.items()}
        ids = sorted(labels, key=kd.get)
        pairs = []
        for i, j in itertools.combinations(ids, 2):
            if kd[i] == kd[j]:
                continue
            pairs.append(
                {
                    "arm": arm,
                    "receptor": rec,
                    "strong": i,
                    "weak": j,
                    "ratio": kd[j] / kd[i],
                    "credit": 1.0
                    if scores[i] < scores[j]
                    else 0.5
                    if scores[i] == scores[j]
                    else 0.0,
                    "matched": all(
                        labels[i][k] == labels[j][k] for k in ("heavy", "formal_charge")
                    ),
                }
            )
        pair_details.extend(pairs)
        for limit, matched in ((1, False), (3, False), (10, False), (3, True)):
            selected = [
                p
                for p in pairs
                if kd[p["weak"]] >= limit * kd[p["strong"]]
                and (not matched or p["matched"])
            ]
            exp = [-np.log10(float(kd[c])) for c in ids]
            pred = [-float(scores[c]) for c in ids]
            summary.append(
                {
                    "arm": arm,
                    "receptor": rec,
                    "min_ratio": limit,
                    "match_charge_size": matched,
                    "pairs": len(selected),
                    "credits": sum(p["credit"] for p in selected),
                    "spearman_all": spearmanr(exp, pred).statistic,
                    "kendall_all": kendalltau(exp, pred).statistic,
                }
            )
    write(path / "ranks.csv", summary)
    write(path / "pairs.csv", pair_details)
    # Fixed first coordinate and self-reference in all three Pro binaries.
    fixed = []
    with tempfile.TemporaryDirectory() as temp:
        pose = Path(temp) / "pose.pdbqt"
        for _, rec, c, protein, box, ligand in entries:
            key = f"Gal3_{rec}_{c}_391781_parent_3000"
            block = (path / "raw" / (key + ".pdbqt")).read_text().split("ENDMDL")[
                0
            ] + "ENDMDL\n"
            pose.write_text(block)
            terms = []
            # A pose printed to 0.001 A can round over a nondecimal box face.
            # Record rejections at the original box. In that case repeat just
            # the score-only parity check at identical coordinates in all arms.
            expanded = False
            for box_attempt in range(2):
                checking = box
                if box_attempt:
                    checking = Path(temp) / "expanded.box"
                    checking.write_text(
                        re.sub(
                            r"(size_[xyz]\s*=\s*)([.\d]+)",
                            lambda m: m[1] + str(float(m[2]) + 0.75),
                            box.read_text(),
                        )
                    )
                    expanded = True
                terms = []
                rejected = []
                for arm in arms[:3]:
                    result = subprocess.run(
                        [
                            getattr(a, arm),
                            "--receptor",
                            protein,
                            "--config",
                            checking,
                            "--ligand",
                            ligand,
                            "--pose_coordinates",
                            pose,
                            "--score_only",
                            "--cpu",
                            "1",
                            "--seed",
                            "391781",
                        ],
                        text=True,
                        capture_output=True,
                        check=False,
                    )
                    if (
                        result.returncode
                        and "outside the grid box" not in result.stderr
                    ):
                        raise ValueError(result.stderr)
                    rejected.append(result.returncode != 0)
                    terms.append(
                        re.findall(r":\s*([-+.\deE]+)\s+\(kcal/mol\)", result.stdout)
                    )
                if not any(rejected):
                    break
                if not all(rejected):
                    raise ValueError(key + " unequal box acceptance")
            if len(terms[0]) != 10 or not terms[0] == terms[1] == terms[2]:
                raise RuntimeError(key)
            fixed.append(
                {
                    "key": key,
                    "sha256": hashlib.sha256(pose.read_bytes()).hexdigest(),
                    "expanded_score_box": expanded,
                    "terms": " ".join(terms[0]),
                }
            )
    write(path / "fixed.csv", fixed)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for x in ("work", "upstream", "upstream-counts", "parent", "always", "adaptive"):
        p.add_argument("--" + x, type=Path, required=True)
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    p.add_argument(
        "--replay-scores",
        action="store_true",
        help="Recompute tables and fixed scores from an existing runs.csv",
    )
    main(p.parse_args())
