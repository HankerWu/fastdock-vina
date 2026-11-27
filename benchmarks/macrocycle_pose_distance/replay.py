"""Locate the first unequal selection on the same observed native minima."""
import argparse
import csv
import itertools
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from benchmarks.concerted_ring.run import command, counts, sha, timed, write


def read_events(log):
    events = []
    for line in log.read_text().splitlines():
        if line.startswith("POSE_EVENT "):
            _, step, energy, *atoms = line.split()
            xyz = np.asarray(atoms, dtype=float).reshape(-1, 4)
            events.append((int(step), float(energy), xyz))
        elif line.startswith("RESIDENTS "):
            step, energy, xyz = events.pop()
            events.append((step, energy, xyz, list(map(float, line.split()[1:]))))
    return events


def distance(a, b, include_markers):
    assert np.array_equal(a[:, 0], b[:, 0])
    diff = a[:, 1:] - b[:, 1:]
    if not include_markers:
        diff = diff[a[:, 0] == 0]
    return np.sqrt(np.sum(diff**2, axis=1).mean())


def update(pool, step, points, include_markers, radius, limit):
    _, en, xyz, _ = points[step]
    near = [s for s in pool if distance(xyz, points[s][2], include_markers) < radius]
    if any(points[s][1] <= en for s in near):
        return pool
    return sorted([s for s in pool if s not in near] + [step], key=lambda s:points[s][1])[:limit]


def compare_stream(traces, key):
    events = {arm: read_events(p) for arm, p in traces.items()}
    points = {arm:{v[0]:v for v in ev} for arm,ev in events.items()}
    pools = {"parent":[], "fixed":[]}
    first = None
    shared = 0
    prefix = True

    def inspect(stage, step):
        nonlocal first
        if first is not None:
            return
        if pools["parent"] != pools["fixed"]:
            # First unequal membership is meaningful only on a common candidate stream.
            assert prefix
            competing = set(previous["parent"]) | {step}
            metric_pairs = []
            for i,j in itertools.combinations(sorted(competing),2):
                a,b = points["parent"][i][2], points["parent"][j][2]
                d0,d1 = distance(a,b,True),distance(a,b,False)
                if (d0 < radius) != (d1 < radius):
                    metric_pairs.append((i,j,d0,d1))
            assert metric_pairs, (stage, step, previous, pools)
            i,j,d0,d1 = metric_pairs[0]
            first = dict(stage=stage,trigger=step,pair_i=i,pair_j=j,old_A=d0,new_A=d1,
                         radius_A=radius,old_residents=" ".join(map(str,pools["parent"])),
                         new_residents=" ".join(map(str,pools["fixed"])))

    radius = 1.0
    for index in range(max(map(len,events.values()))):
        ev = {arm: ls[index] for arm,ls in events.items() if index < len(ls)}
        if len(ev) != 2 or ev["parent"][0:2] != ev["fixed"][0:2] or not np.array_equal(ev["parent"][2], ev["fixed"][2]):
            prefix = False
        if prefix:
            shared += 1
        previous = dict(pools)
        for arm,(step,energy,_,residents) in ev.items():
            pools[arm] = update(pools[arm],step,points[arm],arm=="parent",radius,5)
            assert [points[arm][s][1] for s in pools[arm]] == residents
        if prefix:
            inspect("trajectory",ev["parent"][0])
    # Single trajectory, native cross-trajectory merge has radius 2.
    if first is None:
        assert prefix
        stream = pools["parent"]
        pools = {"parent":[], "fixed":[]}
        radius = 2.0
        for step in stream:
            previous = dict(pools)
            for arm in pools:
                pools[arm] = update(pools[arm],step,points[arm],arm=="parent",radius,5)
            inspect("merge",step)
    return dict(key=key, parent_events=len(events["parent"]),fixed_events=len(events["fixed"]),exact_shared_prefix=shared,
                **(first or dict.fromkeys(("stage","trigger","pair_i","pair_j","old_A","new_A","radius_A","old_residents","new_residents"),"")))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs","production","output","parent","fixed"):
        parser.add_argument(name,type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    with (args.production/"runs.csv").open() as f:
        rows = list(csv.DictReader(f))
    jobs = sorted({(r["pdb"],int(r["seed"]),int(r["cap"])) for r in rows})
    result = []
    for pdb,seed,cap in jobs:
        folder=args.inputs/pdb
        case=(pdb,"",folder/"receptor.pdbqt",folder/"ligand.pdbqt",args.production/(pdb+".box"))
        labels={}
        for arm in ("parent","fixed"):
            key=f"{pdb}_{seed}_{cap}_{arm}"
            log=args.output/(key+".log")
            pose=args.output/(key+".pdbqt")
            timed(command(case,seed,cap,getattr(args,arm).resolve(),pose,True),log)
            prod=args.production/"raw"
            assert sha(pose)==sha(prod/pose.name) and counts(log.read_text())==counts((prod/log.name).read_text())
            pose.unlink()
            labels[arm]=log
        result.append(compare_stream(labels, f"{pdb}_{seed}_{cap}"))
    write(args.output/"summary.csv",result)
