"""Check the ligand coordinates at the first unequal E/gradient evaluation."""
import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

from replay import read_evals, rows


def report(work, data, output):
    dry = (work/"prepared/2X6E.pdbqt").read_bytes()
    full = (work/"template1OL7/trans_model/released/receptor.pdbqt").read_bytes()
    assert full.startswith(dry)
    added = np.array([[float(line[i:i+8]) for i in (30, 38, 46)]
                      for line in full[len(dry):].decode().splitlines()
                      if line.split()[-1] not in ("H","HD")])
    assert len(added) == 101
    center = np.array(json.loads((work/"prepared/2X6E.json").read_text())["center"])
    spacing, nvoxels = 0.375, 70  # ceil(26 / 0.375), even
    lo = center - spacing*nvoxels/2
    hi = center + spacing*nvoxels/2
    offsets = np.array(list(itertools.product((0, 1), repeat=3)))
    results = []
    for row in csv.DictReader((data/"divergences.csv").open()):
        drug, seed, tr = row["drug"], int(row["seed"]), int(row["trajectory"])
        for tag, step, stage in (
            ("first_path", row["first_step"], row["probe_path"]),
            ("first_changed_return", row["first_changed_local_geometry_step"], row["probe_coordinates"]),
        ):
            if not step:
                continue
            logs = [data/"raw"/f"{drug}_{seed}_{arm}_focus.trace" for arm in ("dry", "released")]
            (left,_), (right,_) = (read_evals(p, tr, int(step), stage) for p in logs)
            a, b = next((a,b) for a,b in zip(left,right) if a != b)
            assert a[0] == b[0] and a[2] == b[2] and a[-1] == b[-1]
            xyz = np.array(a[-1]).reshape(-1, 3)
            query = np.clip(xyz, lo, hi)
            index = np.clip(np.floor((query-lo)/spacing).astype(int), 0, nvoxels-1)
            corners = lo + (index[:,None,:] + offsets[None,:,:])*spacing
            # For a clamped boundary coordinate only face/edge nodes may have
            # nonzero value weight. Report conservative minimum over the cell.
            results.append(dict(
                drug=drug, seed=seed, trajectory=tr, probe=tag, step=step, stage=stage, local_eval=a[0],
                n_heavy=len(xyz), n_outside=int(np.count_nonzero((xyz != query).any(axis=1))),
                max_axis_excess_A=np.max(np.abs(xyz-query)),
                nearest_added_direct_A=cdist(xyz, added).min(),
                nearest_added_clamped_A=cdist(query, added).min(),
                nearest_cell_corner_to_added_A=cdist(corners.reshape(-1,3), added).min(),
                dry_energy=a[1], released_energy=b[1], dry_initial=left[0][1], released_initial=right[0][1],
                both_rejected_by_initial_energy_bound=(a[0]>1 and a[1]>=left[0][1] and b[1]>=right[0][1]),
            ))
    rows(output, results)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("data", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    report(a.work, a.data, a.output)
