"""Replay both randomized mechanisms and locate the first different actual rigid proposal."""

import argparse
import os
import sys
from pathlib import Path

import numpy as np

import study

sys.path.insert(0, str(study.HERE.parent / "basin_trace"))
import replay_basin as observer


def first_split(tr, trajectory, step, refs):
    rows = tr[tr["traj"] == trajectory]
    earlier = rows[rows["step"] < step]
    accepted = earlier[earlier["event"] == 5]
    if len(accepted):
        last = accepted[-1]["step"]
        successful = earlier[(earlier["step"] == last) & np.isin(earlier["event"], [5, 6])]
        parent = successful[-1]["xyz"]
    else:
        parent = earlier[earlier["event"] == 1][0]["xyz"]
    branch = rows[rows["step"] == step]
    start = branch[(branch["phase"] == 1) & (branch["local"] == 1) & (branch["event"] == 2)][0]
    returned = branch[branch["event"] == 4][0]
    distance = lambda xyz: min(np.sqrt(np.mean(np.sum((xyz - ref)**2, axis=1))) for ref in refs)
    return parent, start["xyz"], returned["xyz"], dict(
        parent_A=distance(parent), start_A=distance(start["xyz"]), return_A=distance(returned["xyz"]),
        energy_hunt_start=float(start["energy"]), energy_hunt_return=float(returned["energy"]),
        hunt_evals=int(returned["local"]), accepted=int(5 in branch["event"]),
        entered_save=int(6 in branch["event"]), retained_then=int(7 in branch["event"]),
    )


def replay(args):
    out = args.out
    out.mkdir(parents=True, exist_ok=False)
    inputs = study.scan.config(args.aurora, args.brd9, out)
    stages, splits, costs, xyz = [], [], [], {}
    for name, cell in inputs.items():
        for seed in study.SEEDS:
            traces = {}
            refs = None
            for arm in ("pivot", "sham"):
                stem = f"{name}_{seed}_{arm}"
                pool = out / (stem + ".pdbqt")
                trace = out / (stem + ".bin")
                cmd = study.scan.command(args.build / arm / "vina", cell, seed, pool)
                os.environ["FD_BASIN_TRACE"] = str(trace)
                cost = study.scan.replay.invoke(cmd, pool.with_suffix(".log"))
                assert (args.primary / "raw" / pool.name).read_bytes() == pool.read_bytes()
                counts = study.totals(pool.with_suffix(".log"))
                assert counts == study.totals(args.primary / "raw" / pool.with_suffix(".log").name)
                costs.append(dict(drug=name, seed=seed, arm=arm, **cost, **counts,
                                  raw_bytes=trace.stat().st_size, raw_sha256=study.scan.replay.digest(trace)))
                ids, tr = observer.read_trace(trace)
                refs = cell["refs"][:, [list(cell["serials"]).index(j) for j in ids]]
                xyz[name + "_serials"] = ids
                distances = observer.distances(tr, refs)
                labels = observer.stage_labels(tr)
                for stage in np.unique(labels):
                    ids_stage = np.flatnonzero(labels == stage)
                    idx = ids_stage[distances[ids_stage].argmin()]
                    row = tr[idx]
                    xyz[stem + "_" + stage] = row["xyz"].copy()
                    stages.append(dict(drug=name, seed=seed, arm=arm, stage=stage, count=len(ids_stage),
                                       nearest_A=distances[idx], le2=int((distances[ids_stage] <= 2).sum()),
                                       le3=int((distances[ids_stage] <= 3).sum()), trajectory=row["traj"],
                                       step=row["step"], local_eval=row["local"], energy=row["energy"]))
                # Only parent/start/hunt/accepted/saved events are needed to
                # verify a shared prefix; other gradients stream through disk.
                keep = np.isin(tr["event"], (1, 4, 5, 6, 7)) | ((tr["event"] == 2) & (tr["phase"] == 1) & (tr["local"] == 1))
                traces[arm] = tr[keep].copy()
                del tr
                if not args.keep_traces:
                    trace.unlink()
            for trajectory in (1, 2):
                def starts(arm):
                    tr = traces[arm]
                    return tr[(tr["traj"] == trajectory) & (tr["event"] == 2)]
                a, b = starts("sham"), starts("pivot")
                n = min(len(a), len(b))
                idx = np.flatnonzero(np.any(a["xyz"][:n] != b["xyz"][:n], axis=(1, 2)))[0]
                step = a[idx]["step"]
                assert step == b[idx]["step"]
                common = {}
                for arm in ("sham", "pivot"):
                    tr = traces[arm]
                    before = tr[(tr["traj"] == trajectory) & (tr["step"] < step)]
                    common[arm] = before.tobytes()
                    parent, start, returned, diag = first_split(tr, trajectory, step, refs)
                    label = f"{name}_{seed}_{trajectory}_{arm}"
                    xyz[label + "_parent"] = parent
                    xyz[label + "_proposal"] = start
                    xyz[label + "_return"] = returned
                    dist = np.linalg.norm(start - parent, axis=1)
                    diag.update(rms_step_A=np.sqrt(np.mean(dist**2)), fixed_atom_min_A=min(dist),
                                moved_centroid_A=np.linalg.norm((start-parent).mean(axis=0)),
                                max_internal_delta_A=np.max(np.abs(
                                    np.linalg.norm(start[:,None]-start,axis=2)
                                    - np.linalg.norm(parent[:,None]-parent,axis=2))))
                    splits.append(dict(drug=name, seed=seed, trajectory=trajectory, arm=arm, step=step, **diag))
                assert common["pivot"] == common["sham"]
                assert splits[-1]["fixed_atom_min_A"] < 1e-10
                assert abs(splits[-1]["rms_step_A"] - splits[-2]["rms_step_A"]) < 1e-10
    for file, rows in (("stages.csv", stages), ("splits.csv", splits), ("observation_cost.csv", costs)):
        study.scan.replay.geom.write_csv(out, file, rows)
    np.savez_compressed(out / "witnesses.npz", **xyz)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("aurora", "brd9", "primary", "build", "out"):
        p.add_argument(k, type=lambda s: Path(s).resolve())
    p.add_argument("--keep-traces", action="store_true", help="Retain all raw binary coordinate events")
    replay(p.parse_args())
