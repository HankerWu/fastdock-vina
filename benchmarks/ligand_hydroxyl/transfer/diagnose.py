"""Post hoc contact-deletion tests on both poses of every lost Top-1 case."""

import argparse
import copy
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from run_transfer import BENCH

sys.path.insert(0, str(BENCH / "ligand_hydroxyl"))
from model import ACCEPTORS, Receptor, atoms, magnitude, score, topology


def diagnose(work, result):
    cells = pd.read_csv(result / "cells.csv")
    poses = pd.read_csv(result / "poses.csv")
    rows = []
    for (q, r, cap, seed), chunk in cells[cells.engine == "pro"].groupby(
        ["query", "receiver", "cap", "seed"]
    ):
        ranks = chunk.set_index("arm")
        if ranks.loc["native", "top1"] != 1 or ranks.loc["bonded", "top1"] != 0:
            continue
        subset = poses[
            (poses.engine == "pro")
            & (poses["query"] == q)
            & (poses.receiver == r)
            & (poses.cap == cap)
            & (poses.seed == seed)
        ]
        lig = (work / q / "ligand.pdbqt").read_text()
        rec = Receptor((work / r / "receptor.pdbqt").read_text())
        blocks = [
            b
            for b in (work / "runs" / f"{q}_{r}_{cap}_{seed}_pro.pdbqt")
            .read_text()
            .split("ENDMDL")
            if b.strip()
        ]
        _, _, sites = topology(lig)
        for arm in ["native", "bonded"]:
            rank = ranks.loc[arm, "index"]
            pose = subset[subset.model == rank].iloc[0]
            block = blocks[rank - 1]
            xyz = atoms(block)
            nearby = {
                rec.ids[k]
                for o in sites
                for k in rec.tree.query_ball_point(xyz[o].xyz, 3.6)
                if rec.atoms[rec.ids[k]].kind in ACCEPTORS
                and magnitude(xyz[o], rec.atoms[rec.ids[k]]) > 0
            }
            row = dict(
                query=q,
                receiver=r,
                cap=cap,
                seed=seed,
                winner=arm,
                rank=rank,
                rmsd=pose.rmsd,
                native=pose.native,
                bonded=pose.bonded,
                increment=pose.bonded - pose.native,
            )
            if not nearby:
                rows.append(dict(**row, contact="", distance="", remove_delta=0.0))
            for j in sorted(nearby):
                atom = rec.atoms[j]
                masked = copy.copy(rec)
                masked.ids = [i for i in rec.ids if i != j]
                masked.tree = cKDTree([rec.atoms[i].xyz for i in masked.ids])
                delta = score(masked, lig, block)["bonded"] * pose.scale - (
                    pose.bonded - pose.native
                )
                rows.append(
                    dict(
                        **row,
                        contact=atom.residue.strip() + ":" + atom.name,
                        distance=min(
                            np.linalg.norm(xyz[o].xyz - atom.xyz) for o in sites
                        ),
                        remove_delta=delta,
                    )
                )
    pd.DataFrame(rows).to_csv(result / "loss_contacts.csv", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("result", type=Path)
    a = p.parse_args()
    diagnose(a.work, a.result)
