"""Recompute cell QA, paired affinity direction, rank and cost from archived pools."""

import argparse
import itertools
import math
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from statistics import mean

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree
from scipy.stats import kendalltau, spearmanr

from assay import COUNT_KEYS, STATES, read, sha, write


def atoms(block):
    return {
        int(l[6:11]): [float(l[k : k + 8]) for k in (30, 38, 46)]
        for l in block.splitlines()
        if l.startswith(("ATOM  ", "HETATM")) and l.split()[-1] not in ("H", "HD", "HS")
    }


def evaluate(repo, work):
    lookup = {r["compound"]: r for r in read(work / "panel.csv")}
    errors = {
        r["compound"]: r for r in read(repo / "benchmarks/galectin_affinity/table1.csv")
    }
    rows = [r for p in sorted(work.glob("runs-*.csv")) for r in read(p)]
    blocks = defaultdict(dict)
    geometries = []
    costs = []
    prepare = work / "prepared"
    trees = {
        rec: cKDTree(list(atoms((prepare / (rec + ".pdbqt")).read_text()).values()))
        for states in STATES.values()
        for rec in states
    }
    for row in rows:
        arm = row["arm"]
        panel = f"{arm}-{row['seed']}-{row['cap']}"
        rec, cid = row["receptor"], row["compound"]
        assert cid not in blocks[panel, rec]
        blocks[panel, rec][cid] = Decimal(row["top_score"])
        path = work / "raw" / (row["key"] + ".pdbqt")
        assert sha(path) == row["sha256"]
        models = path.read_text().split("ENDMDL")[:-1]
        assert len(models) == int(row["modes"])
        ref = np.loadtxt(prepare / f"{rec}_{cid}.anchor")
        for rank, model in enumerate(models, 1):
            coords = atoms(model)
            xyz = np.array([coords[int(serial)] for serial in ref[:, 0]])
            rmsd = np.mean(((xyz - ref[:, 1:]) ** 2).sum(axis=1)) ** 0.5
            dists = trees[rec].query(list(coords.values()))[0]
            affinity = next(
                l for l in model.splitlines() if l.startswith("REMARK VINA RESULT:")
            ).split()[3]
            if rank == 1:
                assert Decimal(affinity) == blocks[panel, rec][cid]
            geometries.append(
                {
                    "panel": panel,
                    "key": row["key"],
                    "rank": rank,
                    "score": affinity,
                    "donor_ring_A": round(rmsd, 6),
                    "success": int(rmsd <= 2),
                    "minimum_heavy_contact_A": round(min(dists), 6),
                    "atoms_below_2A": int((dists < 2).sum()),
                }
            )
    for panel, rec in blocks:
        cells = [
            r
            for r in rows
            if f"{r['arm']}-{r['seed']}-{r['cap']}" == panel and r["receptor"] == rec
        ]
        keys = {r["key"] for r in cells}
        geom = defaultdict(list)
        for g in geometries:
            if g["key"] in keys:
                geom[g["key"]].append(g)
        hits = [
            sum(any(g["success"] for g in case[:k]) for case in geom.values())
            for k in (1, 5, 9)
        ]
        costs.append(
            dict(
                panel=panel,
                receptor=rec,
                cells=len(cells),
                seconds=round(sum(float(r["seconds"]) for r in cells), 6),
                modes=sum(int(r["modes"]) for r in cells),
                top1_ring=hits[0],
                top5_ring=hits[1],
                oracle_ring=hits[2],
                first_with_sub2A_contact=sum(
                    g[0]["atoms_below_2A"] > 0 for g in geom.values()
                ),
                **{
                    col: sum(int(r[col]) for r in cells) if cells[0][col] else ""
                    for col in COUNT_KEYS
                },
            )
        )
    credit = lambda gap: 1 if gap > 0 else 0 if gap < 0 else 0.5
    edges, results = [], []
    cards = sorted(lookup, key=int)
    size = {i: int(lookup[i]["heavy_atoms"]) for i in cards}
    charge = {
        i: Chem.GetFormalCharge(Chem.MolFromSmiles(lookup[i]["smiles"])) for i in cards
    }
    panels = sorted({p for p, r in blocks if not p.startswith("counts-")})
    for rec in trees:
        blocks["heavy_atoms", rec] = {i: Decimal(-size[i]) for i in cards}
    panels.append("heavy_atoms")
    for target, states in STATES.items():
        kd = {i: Decimal(lookup[i][target + "_Kd_uM"]) for i in cards}
        err = {i: Decimal(errors[i][target + "_error_uM"]) for i in cards}
        for panel in panels:
            for state in (*states, "consensus"):
                used = states if state == "consensus" else (state,)
                pred = {i: mean(blocks[panel, r][i] for r in used) for i in cards}
                x, y = ([float(v[i]) for i in cards] for v in (pred, kd))
                rho, tau = spearmanr(x, y).statistic, kendalltau(x, y).statistic
                subsets = defaultdict(list)
                for j, k in itertools.combinations(cards, 2):
                    a, b = sorted((j, k), key=kd.get)
                    if kd[b] < 3 * kd[a]:
                        continue
                    fold = kd[b] / kd[a]
                    gap = pred[b] - pred[a]
                    row = {
                        "target": target,
                        "state": state,
                        "panel": panel,
                        "strong": a,
                        "weak": b,
                        "fold": float(fold),
                        "gap": str(gap),
                        "credit": credit(gap),
                        "equal_HAC_and_charge": int(
                            size[a] == size[b] and charge[a] == charge[b]
                        ),
                        "triazoles": int("17" not in (a, b)),
                        "tenfold": int(fold >= 10),
                        "measured_error_robust": int(
                            kd[b] - err[b] >= 3 * (kd[a] + err[a])
                        ),
                        "experimental_gap_298": 0.001987204258
                        * 298.15
                        * math.log(float(fold)),
                    }
                    edges.append(row)
                    for name in (
                        "all",
                        "equal_HAC_and_charge",
                        "triazoles",
                        "tenfold",
                        "measured_error_robust",
                    ):
                        if name == "all" or row[name]:
                            subsets[name].append(row)
                for name, sub in subsets.items():
                    results.append(
                        {
                            "target": target,
                            "state": state,
                            "panel": panel,
                            "subset": name,
                            "pairs": len(sub),
                            "correct": sum(r["credit"] for r in sub),
                            "spearman": rho if name == "all" else "",
                            "kendall_tau_b": tau if name == "all" else "",
                            "relative_gap_MAE_298": mean(
                                abs(float(r["gap"]) - r["experimental_gap_298"])
                                for r in sub
                            )
                            if panel != "heavy_atoms"
                            else "",
                        }
                    )
    write(work / "geometry.csv", geometries)
    write(work / "cost.csv", costs)
    write(work / "pairs.csv", edges)
    write(work / "rankings.csv", results)
    # Whole compounds, keeping identical sampled identities across the two arms and receptor states.
    draw = np.random.default_rng(912751).multinomial(
        len(cards), [1 / len(cards)] * len(cards), 4000
    )
    summaries = []
    for candidate_panel in panels:
        if not candidate_panel.startswith("candidate-"):
            continue
        _, seed, cap = candidate_panel.split("-")
        for target, states in STATES.items():
            for state in (*states, "consensus"):
                for name in (
                    "all",
                    "tenfold",
                    "equal_HAC_and_charge",
                    "triazoles",
                    "measured_error_robust",
                ):
                    matched = {}
                    for arm in ("parent", "candidate"):
                        matched[arm] = {
                            (r["strong"], r["weak"]): r["credit"]
                            for r in edges
                            if r["state"] == state
                            and r["target"] == target
                            and r["panel"] == f"{arm}-{seed}-{cap}"
                            and (name == "all" or r[name])
                        }
                    if not matched["parent"] or not matched["candidate"]:
                        continue
                    assert matched["candidate"].keys() == matched["parent"].keys()
                    eligible = list(matched["parent"])
                    weights = np.stack(
                        [
                            draw[:, cards.index(a)] * draw[:, cards.index(b)]
                            for a, b in eligible
                        ],
                        axis=1,
                    )
                    gaps = np.array(
                        [
                            matched["candidate"][p] - matched["parent"][p]
                            for p in eligible
                        ]
                    )
                    n = weights.sum(axis=1)
                    good = n > 0
                    diff = (weights[good] @ gaps) / n[good]
                    lo, hi = np.quantile(diff, (0.025, 0.975))
                    summaries.append(
                        {
                            "target": target,
                            "state": state,
                            "seed": seed,
                            "cap": cap,
                            "subset": name,
                            "pairs": len(eligible),
                            "gains": sum(gaps > 0),
                            "losses": sum(gaps < 0),
                            "estimate_pp": round(float(gaps.mean()) * 100, 6),
                            "lower_pp": round(lo * 100, 6),
                            "upper_pp": round(hi * 100, 6),
                            "valid_draws": good.sum(),
                        }
                    )
    if summaries:
        write(work / "paired_bootstrap.csv", summaries)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("repo", type=Path)
    p.add_argument("work", type=Path)
    a = p.parse_args()
    evaluate(a.repo, a.work)
