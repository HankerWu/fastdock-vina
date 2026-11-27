"""Rank a fixed lattice of rigid-body hypotheses after one optional local step."""

import sys
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
sys.path.insert(0, str(HERE.parent / "rigid_translation_scan"))
from rigid_lookahead import rigid_trial
from scan import pose_batches, select_starts


def translate_search(
    fields, origin, spacing, center, ligand, seed, arm, observe=None, budget=9
):
    if budget not in (9, 24, 48):
        raise ValueError("budget must be 9, 24 or 48")
    heavy, types = ligand[3:]
    proposals = []
    count = 0
    for index, rotated, positions in pose_batches(center, ligand, seed, "lattice"):
        rng = np.random.default_rng(np.random.SeedSequence([seed, index]))
        result = rigid_trial(
            fields, origin, spacing, rotated[heavy], types, positions, arm, rng
        )
        if observe is not None:
            observe(index, rotated, positions, result)

        def candidate(j):
            coords = (
                Rotation.from_rotvec(result["rotvec"][j]).apply(rotated)
                + positions[j]
                + result["translation"][j]
            )
            return float(result["energy"][j]), index, coords

        if budget == 48:
            # Each orientation keeps the raw parent winner and a distinct second start.
            j = int(np.argmin(result["parent_energy"]))
            parent = rotated + positions[j]
            proposals.append((float(result["parent_energy"][j]), index, parent))
            for k in np.argsort(result["energy"], kind="stable"):
                child = candidate(k)
                if not np.allclose(child[2][heavy], parent[heavy], atol=1e-10, rtol=0):
                    proposals.append(child)
                    break
        else:
            proposals.append(candidate(int(np.argmin(result["energy"]))))
        count += len(positions)
    starts = (
        select_starts(proposals, heavy)
        if budget == 9
        else sorted(proposals, key=lambda item: (item[0], item[1]))
    )
    return proposals, starts, count
