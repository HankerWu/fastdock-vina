"""Replay original full assay panels on the exact same saved native first poses."""

import argparse
from decimal import Decimal
from functools import lru_cache
from itertools import combinations
import hashlib
import json
from math import log10
from pathlib import Path
import sys
import time

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT.parent
sys.path[:0] = [
    str(ROOT),
    str(BENCH / "prediction_metrics"),
    str(BENCH / "terminal_rotor_score"),
]
from model import AmideReceptor
from coupled import CoupledAmides
from affinity import ASSAYS, UNIT_M, evaluate
from assay import WEIGHT, aggregate, assay_inputs, read, save

ARMS = ("native", "ensemble", "coupled", "minimum", "uniform", "flip_all")


def chemical_controls(query):
    line = next(
        l
        for l in query.splitlines()
        if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX")
    )
    mol = Chem.MolFromSmiles(line.split()[2])
    if mol is None:
        raise ValueError("Unreadable prepared ligand chemistry")
    return mol.GetNumHeavyAtoms(), Chem.GetFormalCharge(mol)


def run(out):
    out.mkdir(parents=True, exist_ok=True)
    src = {
        (r["target"], r["ligand"], r["receptor"], int(r["seed"])): r
        for r in read(BENCH / "terminal_rotor_score/results/poses.csv")
    }
    rec_cost, poses, sites = [], [], []

    @lru_cache(None)
    def receptor(blob):
        start = time.perf_counter()
        rec = AmideReceptor(blob.decode())
        mid = time.perf_counter()
        graph = CoupledAmides(rec)
        rec_cost.append(
            dict(
                sha256=hashlib.sha256(blob).hexdigest(),
                sites=len(rec.sites),
                edges=len(graph.edges),
                components=len(graph.components),
                max_width=max(c.width for c in graph.components),
                independent_seconds=mid - start,
                graph_seconds=time.perf_counter() - mid,
            )
        )
        return rec, graph

    ligand_cost = joint_cost = 0.0
    for target, ligand, measured, receive, seed, protein, query, pool in assay_inputs():
        row = src.pop((target, ligand, receive, seed))
        pose = pool.decode().split("ENDMDL", 1)[0]
        native, scale = (
            float(row["fixed"]),
            1 / (1 + WEIGHT * float(row["fixed_rotors"])),
        )
        if abs(native - scale * float(row["saved_inter"])) > 1.0e-10 or float(
            pose.split("REMARK VINA RESULT:")[1].split()[0]
        ) != float(row["saved_score"]):
            raise ValueError("Archived first pose or normalization mismatch")
        rec, graph = receptor(protein)
        text = query.decode()
        nheavy, charge = chemical_controls(text)
        start = time.perf_counter()
        terms = rec.score(text, pose, detail=True)
        mid = time.perf_counter()
        terms["coupled"] = graph.correction(terms["sites"])
        joint_cost += time.perf_counter() - mid
        ligand_cost += mid - start
        key = dict(target=target, ligand=ligand, receptor=receive, seed=seed)
        poses.append(
            dict(
                **key,
                measurement=measured,
                heavy_atoms=nheavy,
                formal_charge=charge,
                scale=scale,
                native=native,
                **{arm: native + scale * terms[arm] for arm in ARMS[1:]},
                active_sites=terms["active_sites"],
            )
        )
        for site, apo, delta, corr in terms["sites"]:
            sites.append(
                dict(
                    **key,
                    site=site.strip(),
                    apo_delta=apo,
                    ligand_delta=delta,
                    independent_correction=corr,
                )
            )
    if src or len(poses) != 234:
        raise ValueError("Incomplete paired assay panel")
    save(out / "poses.csv", poses)
    save(out / "contacts.csv", sites)
    save(out / "receptors.csv", rec_cost)
    (out / "cost.json").write_text(
        json.dumps(
            dict(
                poses=len(poses),
                score_seconds=ligand_cost,
                coupled_extra_seconds=joint_cost,
            ),
            indent=2,
        )
        + "\n"
    )

    ligands, stats, paired, edges = [], [], [], []
    for target in sorted({r["target"] for r in poses}):
        group = [r for r in poses if r["target"] == target]
        assay, endpoint, unit = ASSAYS[target]
        for receptor_pool in ["mean", *sorted({r["receptor"] for r in group})]:
            rows = [
                r
                for r in group
                if receptor_pool == "mean" or r["receptor"] == receptor_pool
            ]
            values = {arm: aggregate(rows, arm) for arm in ARMS}
            panel = []
            for name in sorted(values["native"]):
                runs = [r for r in rows if r["ligand"] == name]
                props = {
                    (r["measurement"], r["heavy_atoms"], r["formal_charge"])
                    for r in runs
                }
                if len(props) != 1:
                    raise ValueError(
                        f"Chemistry or measurement differs across prepared {target}/{name}"
                    )
                measured, heavy, charge = props.pop()
                record = dict(
                    target=target,
                    assay=assay,
                    budget="archived",
                    pose_rule="fixed_native_first",
                    endpoint=endpoint,
                    unit=unit,
                    K=Decimal(str(measured)),
                    pK=-log10(measured * UNIT_M[unit]),
                    receptor=receptor_pool,
                    ligand=name,
                    heavy_atoms=heavy,
                    formal_charge=charge,
                    runs=len(runs),
                    **{k: v[name] for k, v in values.items()},
                )
                ligands.append(record)
                panel.append(record)
            for arm in ARMS[1:]:
                # The same draw matrix is used for all comparisons on this panel.
                s, diff = evaluate(
                    panel, np.random.default_rng(43103), 10000, arms=("native", arm)
                )
                stats.extend(
                    dict(receptor=receptor_pool, **r)
                    for r in (s if arm == "ensemble" else s[1:])
                )
                paired.append(dict(receptor=receptor_pool, candidate=arm, **diff))
            for x, y in combinations(panel, 2):
                if x["K"] == y["K"]:
                    continue
                strong, weak = sorted((x, y), key=lambda r: r["K"])
                edges.append(
                    dict(
                        target=target,
                        receptor=receptor_pool,
                        strong=strong["ligand"],
                        weak=weak["ligand"],
                        fold=float(weak["K"] / strong["K"]),
                        matched=int(
                            x["formal_charge"] == y["formal_charge"]
                            and x["heavy_atoms"] == y["heavy_atoms"]
                        ),
                        **{
                            arm: (np.sign(weak[arm] - strong[arm]) + 1) / 2
                            for arm in ARMS
                        },
                    )
                )
    save(out / "ligands.csv", ligands)
    save(out / "ranks.csv", stats)
    save(out / "differences.csv", paired)
    save(out / "pair_credit.csv", edges)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    run(p.parse_args().output)
