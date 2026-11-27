"""Coupled-state rescoring and exact score-gap bounds for any common apo-state prior."""

import argparse
from collections import defaultdict
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
import pandas as pd
from coupled import CoupledAmides
from evaluate import KEY
from model import AmideReceptor


def analyze(work, report):
    poses = pd.read_csv(report / "poses.csv")
    changes = pd.read_csv(report / "site_changes.csv")
    details = defaultdict(list)
    for row in changes.to_dict("records"):
        key = tuple(row[k] for k in KEY) + (row["model"],)
        details[key].append(
            tuple(
                row[k] for k in ("residue", "apo_delta", "ligand_delta", "correction")
            )
        )
    arms, cells, timings = [], [], []
    for receiver, panel in poses.groupby("receiver"):
        wall, cpu = perf_counter(), process_time()
        system = CoupledAmides(
            AmideReceptor((work / receiver / "receptor.pdbqt").read_text())
        )
        setup = perf_counter() - wall
        prep_cpu = process_time() - cpu
        for pool, frame in panel.groupby(KEY):
            info = dict(zip(KEY, pool))
            rows = list(frame.sort_values("model").itertuples())
            native = rows[0]
            baseline = np.zeros(len(system.index))
            candidates = []
            for row in rows:
                sites = details[pool + (row.model,)]
                field = np.zeros(len(system.index))
                for residue, _, dl, _ in sites:
                    field[system.index[residue]] = dl
                if row.model == 1:
                    baseline = field.copy()
                # A common prior gives a log-sum-exp in which the ratio is a normalized
                # weighted average of exp(-[field_j-field_1].s/RT), s_i in {0,1}.
                # This bound also holds when the receptor states are statistically coupled.
                gap = row.native - native.native
                delta = field - baseline
                assert row.scale == native.scale
                lower = gap + row.scale * np.minimum(0.0, delta).sum()
                upper = gap + row.scale * np.maximum(0.0, delta).sum()
                joint = row.native + row.scale * system.correction(sites)
                candidates.append(
                    dict(
                        **info,
                        model=row.model,
                        rmsd_A=row.rmsd_A,
                        native=row.native,
                        factorized=row.ensemble,
                        joint=joint,
                        native_gap=gap,
                        min_possible_gap=lower,
                        max_possible_gap=upper,
                    )
                )
            top = candidates[0]
            assert top["model"] == 1
            for x in candidates:
                for arm in ("factorized", "joint"):
                    assert (
                        x["min_possible_gap"] - 1.0e-10
                        <= x[arm] - top[arm]
                        <= x["max_possible_gap"] + 1.0e-10
                    )
            challengers = candidates[1:]
            nearest = min((r["min_possible_gap"] for r in challengers), default=np.inf)
            ranked = sorted(candidates, key=lambda r: (r["joint"], r["model"]))
            fact = sorted(candidates, key=lambda r: (r["factorized"], r["model"]))
            first = ranked[0]
            best = min(x["rmsd_A"] for x in ranked)
            rms = [x["rmsd_A"] for x in ranked]
            possible = [x for x in challengers if x["min_possible_gap"] <= 0.001]
            cells.append(
                dict(
                    **info,
                    count=len(ranked),
                    any_changed_order=[x["model"] for x in ranked]
                    != [x["model"] for x in fact],
                    selected=first["model"],
                    native_top1_A=top["rmsd_A"],
                    top1_A=rms[0],
                    best5_A=min(rms[:5]),
                    oracle_A=best,
                    regret_A=rms[0] - best,
                    min_challenger_lower_gap=nearest,
                    guaranteed_top1=not possible,
                    geometrically_better_possible=any(
                        x["rmsd_A"] < top["rmsd_A"] for x in possible
                    ),
                    optimistically_best_A=min(
                        [x["rmsd_A"] for x in possible] + [top["rmsd_A"]]
                    ),
                    unconstrained_challengers=len(possible),
                )
            )
            arms.extend(candidates)
        timings.append(
            {
                "receiver": receiver,
                "baseline_and_joint_setup_s": setup,
                "baseline_and_joint_setup_cpu_s": prep_cpu,
                "score_and_bounds_s": perf_counter() - wall - setup,
                "score_and_bounds_cpu_s": process_time() - cpu - prep_cpu,
                "components": len(system.components),
                "max_width": max((c.width for c in system.components), default=0),
                "sites": len(system.index),
                "nonzero_couplings": len(system.edges),
            }
        )
    for name, data in [
        ("joint_poses", arms),
        ("joint_cells", cells),
        ("joint_cost", timings),
    ]:
        pd.DataFrame(data).to_csv(report / f"{name}.csv", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("report", type=Path)
    args = p.parse_args()
    analyze(args.work, args.report)
