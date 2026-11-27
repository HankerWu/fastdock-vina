"""Reconstruct finite crystal-assisted rigid witnesses and test their receptor contacts."""
import argparse
import json
from pathlib import Path
import tempfile
import time

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

from decompose import BENCH, PANELS, by, published_modes, read_csv, reference_reader, table, unpack


def fits(x, refs):
    xbar, ybar = x.mean(axis=0), refs.mean(axis=1)
    local = x - xbar
    targets = refs - ybar[:, None]
    u, _, vt = np.linalg.svd(local.T @ targets)
    u[:, :, -1] *= np.linalg.det(u @ vt)[:, None]
    rotation = u @ vt
    translation = ybar - xbar @ rotation
    pos = x @ rotation + translation[:, None]
    # The correspondence that generated a placement need not minimize its frame error.
    squared = np.square(pos[:, None] - refs[None]).sum(axis=-1).mean(axis=-1)
    return pos, np.sqrt(squared.min(axis=1)), rotation, translation


def receptor_xyz(work, prepared, name):
    lines = (work / prepared[name]["receptor"]).read_text().splitlines()
    return np.array([[float(line[j:j+8]) for j in (30, 38, 46)] for line in lines
                     if line.startswith(("ATOM  ", "HETATM")) and line.split()[-1] not in ("H", "HD", "HS")])


def analyze(panel, work, results):
    parent = BENCH / panel / "pose"
    unpack(parent, work)
    refs, xyz = reference_reader(panel, parent, work)
    previous, files = published_modes(panel, parent, work)
    prepared = json.loads((work/"inputs.json").read_text())["cells"]
    receptor = {name: receptor_xyz(work, prepared, name) for name in refs}
    trees = {name: cKDTree(pos) for name, pos in receptor.items()}
    old_contacts = {(r["stem"], int(r["rank"])): float(r["min_receptor_distance_A"])
                    for r in read_csv(parent/"results/chemistry.csv")}
    if panel == "rigid_pose_coverage":
        for r in read_csv(parent/"results/extra_trajectory_modes.csv"):
            old_contacts[(r["case"]+"_pro-extra", int(r["rank"]))] = float(r["min_receptor_distance_A"])
    rows, details = [], []
    max_contact_error, max_replay_error = 0., 0.
    for stem, records in by(previous, "stem").items():
        drug = records[0]["drug"]
        cell, tree = refs[drug], trees[drug]
        blocks = (work/"raw"/(files[stem]+".pdbqt")).read_text().split("ENDMDL")[:-1]
        witness = None
        count_near, count_valid, count_fits = 0, 0, 0
        for block, r in zip(blocks, records, strict=True):
            rank = int(r["rank"])
            x = xyz(block, cell["serials"])
            original_contact = tree.query(x)[0].min()
            max_contact_error = max(max_contact_error, abs(original_contact-old_contacts[stem, rank]))
            placed, score, rot, shift = fits(x, cell["refs"])
            closest = tree.query(placed)[0].min(axis=1)
            finite = np.flatnonzero(closest >= 2)
            count_fits += len(placed)
            count_near += int(score.min() <= 2)
            count_valid += int(len(finite) > 0 and score[finite].min() <= 2)
            best = int(finite[score[finite].argmin()]) if len(finite) else None
            details.append(dict(stem=stem, rank=rank, rmsd_A=r["rmsd_A"],
                                original_contact_A=original_contact, maps=len(score), no_short_maps=len(finite),
                                min_A=score.min(), best_no_short_A=score[best] if best is not None else "",
                                best_no_short_map=best if best is not None else "",
                                best_no_short_contact_A=closest[best] if best is not None else ""))
            if best is not None and (witness is None or score[best] < witness[0]):
                witness = (score[best], rank, best, closest[best], rot[best], shift[best], x)
        row = dict(stem=stem, drug=drug, seed=records[0]["seed"], arm=records[0]["arm"],
                   released_oracle_A=min(float(r["rmsd_A"]) for r in records), n=len(records), maps=count_fits,
                   modes_le2=count_near, modes_le2_no_short=count_valid)
        if witness is None:
            rows.append(row | dict(best_no_short_A="", rank="", mapping="", min_contact_A="", rotation="", shift=""))
            continue
        rms, rank, mapping, near, rot, shift, x = witness
        restored = x @ rot + shift
        brute = cdist(restored, receptor[drug]).min()
        direct = np.linalg.norm(restored[None]-cell["refs"], axis=(1, 2)).min()/np.sqrt(len(x))
        max_replay_error = max(max_replay_error, abs(brute-near), abs(direct-rms))
        rows.append(row | dict(best_no_short_A=rms, rank=rank, mapping=mapping, min_contact_A=near,
                               rotation=json.dumps(rot.tolist(), separators=(",", ":")),
                               shift=json.dumps(shift.tolist(), separators=(",", ":"))))
    if max(max_contact_error, max_replay_error) > 1e-8:
        raise ValueError("Original pose contact or rigid witness did not reproduce")
    table(results/(panel+"_placements.csv"), rows)
    table(results/(panel+"_placement_modes.csv"), details)
    return dict(panel=panel, pools=len(rows), modes=len(details), fits=sum(r["maps"] for r in rows),
                near_pools=sum(r["modes_le2"] > 0 for r in rows),
                near_no_short_pools=sum(r["modes_le2_no_short"] > 0 for r in rows),
                severe_pools=sum(r["released_oracle_A"] > 5 for r in rows),
                severe_near_no_short=sum(r["released_oracle_A"] > 5 and r["modes_le2_no_short"] > 0 for r in rows),
                max_original_contact_error_A=max_contact_error, max_witness_replay_error_A=max_replay_error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)
    receipts = []
    for panel in PANELS:
        with tempfile.TemporaryDirectory(prefix="rigid_witness_") as temp:
            tic, cpu = time.perf_counter(), time.process_time()
            result = analyze(panel, Path(temp), args.out)
            receipts.append(result | dict(wall_s=time.perf_counter()-tic, cpu_s=time.process_time()-cpu))
    table(args.out/"placement_check.csv", receipts)


if __name__ == "__main__":
    main()
