"""Recompute structural, physical and within-assay endpoints from the stored poses."""

import argparse
import itertools
import json
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from statistics import mean

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import kendalltau, spearmanr

from run import HERE, cases, read, write
from geometry import coord_blocks


def bace_reference(work, out, runs):
    # Reuse the registered MOL2-to-deposition mapping and rigid receptor fit.
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        (tmp / "raw").mkdir()
        for r in runs:
            if r["pdb"] == "5QCT":
                (tmp / "raw" / ("single-" + r["key"] + "-0.pdbqt")).symlink_to(
                    (out / "raw" / (r["key"] + ".pdbqt")).resolve()
                )
        subprocess.run(
            [
                sys.executable,
                HERE.parent / "macrocycle_labels/geometry.py",
                "--upstream-data",
                work / "tutorial",
                "--structure",
                work / "bace/5QCT.pdb",
                "--out",
                tmp,
            ],
            check=True,
        )
        return {
            (r["engine"], int(r["rank"])): float(r["rmsd"])
            for r in read(tmp / "geometry.csv")
        }


def structural(work, out, runs):
    bace = bace_reference(work, out, runs) if any(r["pdb"] == "5QCT" for r in runs) else {}
    modes, summaries = [], []
    for pdb, cid, receptor, ligand, _ in cases(work, "macro"):
        initial = next(coord_blocks(ligand))
        ends = [i for i, (ty, _) in initial.items() if ty == "CG0"]
        dummies = [i for i, (ty, _) in initial.items() if ty == "G0"]
        assert len(ends) == len(dummies) == 2

        def closest(i):
            return min(
                dummies, key=lambda j: np.linalg.norm(initial[i][1] - initial[j][1])
            )

        pairs = [(i, closest(i)) for i in ends]
        assert len(set(j for _, j in pairs)) == 2
        rec = [
            v
            for ty, v in next(coord_blocks(receptor)).values()
            if ty not in ("HD", "H", "HS")
        ]
        tree = cKDTree(rec)
        if pdb != "5QCT":
            d = work / ("hsp" if pdb == "1YET" else "extra") / pdb
            meta = json.loads((d / "preparation.json").read_text())
            references = np.asarray(meta["mapped"])
        for r in (r for r in runs if r["pdb"] == pdb):
            pool = []
            for rank, pose in enumerate(
                coord_blocks(out / "raw" / (r["key"] + ".pdbqt")), 1
            ):
                if pdb == "5QCT":
                    rmsd = bace[r["key"], rank]
                else:
                    xyz = np.array([pose[i][1] for i in meta["serial"]])
                    rmsd = np.sqrt(
                        np.mean(np.sum((xyz - references) ** 2, axis=-1), axis=1)
                    ).min()
                heavy = [
                    v
                    for ty, v in pose.values()
                    if ty not in ("HD", "H", "HS") and not ty.startswith("G")
                ]
                dist = tree.query(heavy)[0]
                item = dict(
                    key=r["key"],
                    rank=rank,
                    rmsd_A=round(float(rmsd), 6),
                    closure_bond_A=round(
                        float(np.linalg.norm(pose[ends[0]][1] - pose[ends[1]][1])), 6
                    ),
                    max_dummy_error_A=round(
                        float(
                            max(
                                np.linalg.norm(pose[i][1] - pose[j][1])
                                for i, j in pairs
                            )
                        ),
                        6,
                    ),
                    closest_contact_A=round(float(min(dist)), 6),
                    ligand_atoms_within_2A=int(sum(dist < 2)),
                )
                modes.append(item)
                pool.append(item)
            assert len(pool) == int(r["modes"])
            rmsds = [p["rmsd_A"] for p in pool]
            summaries.append(
                dict(
                    **r,
                    top1=rmsds[0],
                    best5=min(rmsds[:5]),
                    oracle=min(rmsds),
                    top1_clashes=pool[0]["ligand_atoms_within_2A"],
                )
            )
    write(out / "models.csv", modes)
    write(out / "structures.csv", summaries)


def assay(out, runs):
    panel = read(HERE.parent / "search_backtracking/assay/panel.csv")
    pairs, summary = [], []
    for (pdb, arm), rr in itertools.groupby(
        sorted(runs, key=lambda r: (r["pdb"], r["arm"])),
        key=lambda r: (r["pdb"], r["arm"]),
    ):
        group = list(rr)
        s = {r["compound"]: float(r["first_score"]) for r in group}
        entries = sorted(
            [r for r in panel if r["pdb"] == pdb], key=lambda r: float(r["kd_nm"])
        )
        assert set(s) == {r["cid"] for r in entries}
        edges = []
        for a, b in itertools.combinations(entries, 2):
            ratio = float(b["kd_nm"]) / float(a["kd_nm"])
            if ratio < 3:
                continue
            diff = s[b["cid"]] - s[a["cid"]]
            edges.append(
                dict(
                    pdb=pdb,
                    target=a["target"],
                    arm=arm,
                    strong=a["compound"],
                    weak=b["compound"],
                    ratio=ratio,
                    predicted_gap=round(diff, 6),
                    credit=(diff > 0) + 0.5 * (diff == 0),
                    same_size_charge=(a["heavy"], a["charge"])
                    == (b["heavy"], b["charge"]),
                    same_scaffold=a["scaffold"] == b["scaffold"],
                    abs_error_kcal=abs(diff - 0.00198720426 * 298.15 * math.log(ratio)),
                )
            )
        pairs += edges
        labels = [float(r["kd_nm"]) for r in entries]
        scores = [s[r["cid"]] for r in entries]
        for subset in ("all", "same_size_charge", "same_scaffold"):
            sub = [e for e in edges if subset == "all" or e[subset]]
            summary.append(
                dict(
                    pdb=pdb,
                    target=entries[0]["target"],
                    arm=arm,
                    subset=subset,
                    n=len(sub),
                    correct=sum(e["credit"] for e in sub),
                    gap_mae=mean(e["abs_error_kcal"] for e in sub) if sub else "",
                    all_rho=spearmanr(scores, labels).statistic,
                    all_tau=kendalltau(scores, labels).statistic,
                    joint=sum(int(r["joint"]) for r in group),
                    seconds=sum(int(r["elapsed_ns"]) for r in group) / 1e9,
                )
            )
    write(out / "kd_pairs.csv", pairs)
    write(out / "kd_summary.csv", summary)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("study", choices=("kd", "macro"))
    p.add_argument("work", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    runs = read(a.output / "runs.csv")
    for r in runs:
        path = a.output / "raw" / r["key"]
        scores = re.findall(
            r"^REMARK VINA RESULT:\s+([-+.\d]+)",
            path.with_suffix(".pdbqt").read_text(),
            re.M,
        )
        assert len(scores) == int(r["modes"]) and scores[0] == r["first_score"]
        log = path.with_suffix(
            ".counts" if r["arm"] == "upstream" else ".log"
        ).read_text()
        for k in (
            "hunt_evals",
            "second_evals",
            "explicit_evals",
            "hunt_calls",
            "second_calls",
            "explicit_calls",
        ):
            assert re.search(r"\b" + k + r"=(\d+)", log)[1] == r[k]
    if a.study == "kd":
        assay(a.output, runs)
    else:
        structural(a.work, a.output, runs)


if __name__ == "__main__":
    main()
